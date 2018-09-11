/*
 * Copyright (C) 2015, Simon Fuhrmann, Fabian Langguth
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include <stdexcept>
#include <iostream>

#include "util/timer.h"
#include "math/matrix.h"
#include "math/matrix_tools.h"
#include "sfm/ba_linear_solver.h"
#include "sfm/ba_cholesky.h"
#include "sfm/ba_conjugate_gradient.h"

SFM_NAMESPACE_BEGIN
SFM_BA_NAMESPACE_BEGIN

namespace
{
    /*
     * Inverts a symmetric, positive definite matrix with NxN bocks on its
     * diagonal using Cholesky decomposition. All other entries must be zero.
     */
    void
    invert_block_matrix_NxN_inplace (SparseMatrix<double>* A, int blocksize)
    {
        if (A->num_rows() != A->num_cols())
            throw std::invalid_argument("Block matrix must be square");
        if (A->num_non_zero() != A->num_rows() * blocksize)
            throw std::invalid_argument("Invalid number of non-zeros");

        int const bs2 = blocksize * blocksize;
        std::vector<double> matrix_block(bs2);
        for (double* iter = A->begin(); iter != A->end(); )
        {
            double* iter_backup = iter;
            for (int i = 0; i < bs2; ++i)
                matrix_block[i] = *(iter++);

            cholesky_invert_inplace(matrix_block.data(), blocksize);

            iter = iter_backup;
            for (int i = 0; i < bs2; ++i)
                if (std::isfinite(matrix_block[i]))
                    *(iter++) = matrix_block[i];
                else
                    *(iter++) = 0.0;
        }
    }

    /*
     * Inverts a matrix with 3x3 bocks on its diagonal. All other entries
     * must be zero. Reading blocks is thus very efficient.
     */
    void
    invert_block_matrix_3x3_inplace (SparseMatrix<double>* A)
    {
        if (A->num_rows() != A->num_cols())
            throw std::invalid_argument("Block matrix must be square");
        if (A->num_non_zero() != A->num_rows() * 3)
            throw std::invalid_argument("Invalid number of non-zeros");

        for (double* iter = A->begin(); iter != A->end(); )
        {
            double* iter_backup = iter;
            math::Matrix<double, 3, 3> rot;
            for (int i = 0; i < 9; ++i)
                rot[i] = *(iter++);

            double det = math::matrix_determinant(rot);
            if (MATH_DOUBLE_EQ(det, 0.0))
                continue;

            rot = math::matrix_inverse(rot, det);
            iter = iter_backup;
            for (int i = 0; i < 9; ++i)
                *(iter++) = rot[i];
        }
    }

    /*
     * Computes for a given matrix A the square matrix A^T * A for the
     * case that block columns of A only need to be multiplied with itself.
     * Becase the resulting matrix is symmetric, only about half the work
     * needs to be performed.
     */
    void
    matrix_block_column_multiply (SparseMatrix<double> const& A,
        std::size_t block_size, SparseMatrix<double>* B)
    {
        SparseMatrix<double>::Triplets triplets;
        triplets.reserve(A.num_cols() * block_size);
        for (std::size_t block = 0; block < A.num_cols(); block += block_size)
        {
            std::vector<DenseVector<double>> columns(block_size);
            for (std::size_t col = 0; col < block_size; ++col)
                A.column_nonzeros(block + col, &columns[col]);
            for (std::size_t col = 0; col < block_size; ++col)
            {
                double dot = columns[col].dot(columns[col]);
                triplets.emplace_back(block + col, block + col, dot);
                for (std::size_t row = col + 1; row < block_size; ++row)
                {
                    dot = columns[col].dot(columns[row]);
                    triplets.emplace_back(block + row, block + col, dot);
                    triplets.emplace_back(block + col, block + row, dot);
                }
            }
        }
        B->allocate(A.num_cols(), A.num_cols());
        B->set_from_triplets(triplets);
    }
}

LinearSolver::Status
LinearSolver::solve (SparseMatrixType const& jac_cams,
    SparseMatrixType const& jac_points,
    DenseVectorType const& vector_f,
    DenseVectorType* delta_x)
{
    bool const has_jac_cams = jac_cams.num_rows() > 0;
    bool const has_jac_points = jac_points.num_rows() > 0;

    /* Select solver based on bundle adjustment mode. */
    if (has_jac_cams && has_jac_points)
        return this->solve_schur(jac_cams, jac_points, vector_f, delta_x);
    else if (has_jac_cams && !has_jac_points)
        return this->solve(jac_cams, vector_f, delta_x, 0);
    else if (!has_jac_cams && has_jac_points)
        return this->solve(jac_points, vector_f, delta_x, 3);
    else
        throw std::invalid_argument("No Jacobian given");
}

LinearSolver::Status
LinearSolver::solve_schur (SparseMatrixType const& jac_cams,
    SparseMatrixType const& jac_points,
    DenseVectorType const& values, DenseVectorType* delta_x)
{
    /*
     *   雅阁比矩阵：
     *           J = [Jc Jp]
     *   Jc是与相机相关的模块，Jp是与三维点相关的模块。
     *   正规方程
     *          (J^TJ + lambda*I)delta_x = J^T(x - F)
     *   进一步写为
     *   [ Jcc+ lambda*Icc   Jcp            ][delta_c]= [v]
     *   [ Jxp               Jpp+lambda*Ipp ][delta_p]  [w]
     *
     *   B = Jcc, E = Jcp, C = Jpp
     *  其中 Jcc = Jc^T* Jc, Jcx = Jc^T*Jx, Jxc = Jx^TJc, Jxx = Jx^T*Jx
     *      v = Jc^T(F-x), w = Jx^T(F-x), deta_x = [delta_c; delta_p]
     */

    // 误差向量
    DenseVectorType const& F = values;
    // 关于相机的雅阁比矩阵
    SparseMatrixType const& Jc = jac_cams;
    // 关于三维点的雅阁比矩阵
    SparseMatrixType const& Jp = jac_points;
    SparseMatrixType JcT = Jc.transpose();
    SparseMatrixType JpT = Jp.transpose();

    // 构造正规方程
    SparseMatrixType B, C;
    // B = Jc^T* Jc
    matrix_block_column_multiply(Jc, this->opts.camera_block_dim, &B);
    // C = Jp^T*Jp
    matrix_block_column_multiply(Jp, 3, &C);
    // E = Jc^T*Jp
    SparseMatrixType E = JcT.multiply(Jp);

    /* Assemble two values vectors. */
    DenseVectorType v = JcT.multiply(F);
    DenseVectorType w = JpT.multiply(F);
    v.negate_self();
    w.negate_self();

    /* 以矩阵B和C的对角元素重新构建对角阵*/
//    SparseMatrixType B_diag = B.diagonal_matrix();
//    SparseMatrixType C_diag = C.diagonal_matrix();

    /* 添加信赖域 */
    C.mult_diagonal(1.0 + 1.0 / this->opts.trust_region_radius);
    B.mult_diagonal(1.0 + 1.0 / this->opts.trust_region_radius);

    /* 求解C矩阵的逆C = inv(Jx^T+Jx + lambda*Ixx)*/
    invert_block_matrix_3x3_inplace(&C);

    /* 计算S矩阵的Schur补用于高斯消元. */
    SparseMatrixType ET = E.transpose();

    // S = (Jcc+lambda*Icc) - Jc^T*Jx*inv(Jxx+ lambda*Ixx)*Jx^T*Jc
    SparseMatrixType S = B.subtract(E.multiply(C).multiply(ET));
    // rhs = v -  Jc^T*Jx*inv(Jxx+ lambda*Ixx)*w
    DenseVectorType rhs = v.subtract(E.multiply(C.multiply(w)));

    /* Compute pre-conditioner for linear system. */
    //SparseMatrixType precond = S.diagonal_matrix();
    //precond.cwise_invert();
    SparseMatrixType precond = B;
    invert_block_matrix_NxN_inplace(&precond, this->opts.camera_block_dim);

    /* 用共轭梯度法求解相机参数. */
    DenseVectorType delta_y(Jc.num_cols());
    typedef sfm::ba::ConjugateGradient<double> CGSolver;
    CGSolver::Options cg_opts;
    cg_opts.max_iterations = this->opts.cg_max_iterations;
    cg_opts.tolerance = 1e-20;
    CGSolver solver(cg_opts);
    CGSolver::Status cg_status;
    cg_status = solver.solve(S, rhs, &delta_y, &precond);

    Status status;
    status.num_cg_iterations = cg_status.num_iterations;
    switch (cg_status.info)
    {
        case CGSolver::CG_CONVERGENCE:
            status.success = true;
            break;
        case CGSolver::CG_MAX_ITERATIONS:
            status.success = true;
            break;
        case CGSolver::CG_INVALID_INPUT:
            std::cout << "BA: CG failed (invalid input)" << std::endl;
            status.success = false;
            return status;
        default:
            break;
    }

    /* 将相机参数带入到第二个方程中，求解三维点的参数. */
    /*E= inv(Jp^T Jp) (JpT.multiply(F)-Jc^T * Jp * delta_y)*/
    DenseVectorType delta_z = C.multiply(w.subtract(ET.multiply(delta_y)));

    /* Fill output vector. */
    std::size_t const jac_cam_cols = Jc.num_cols();
    std::size_t const jac_point_cols = Jp.num_cols();
    std::size_t const jac_cols = jac_cam_cols + jac_point_cols;

    if (delta_x->size() != jac_cols)
        delta_x->resize(jac_cols, 0.0);
    for (std::size_t i = 0; i < jac_cam_cols; ++i)
        delta_x->at(i) = delta_y[i];
    for (std::size_t i = 0; i < jac_point_cols; ++i)
        delta_x->at(jac_cam_cols + i) = delta_z[i];

    /* Compute predicted error decrease */
//    status.predicted_error_decrease = 0.0;
//    status.predicted_error_decrease += delta_y.dot(B_diag.multiply(
//        delta_y).multiply(1.0 / this->opts.trust_region_radius).add(v));
//    status.predicted_error_decrease += delta_z.dot(C_diag.multiply(
//        delta_z).multiply(1.0 / this->opts.trust_region_radius).add(w));

    return status;
}

LinearSolver::Status
LinearSolver::solve (SparseMatrixType const& J,
    DenseVectorType const& vector_f,
    DenseVectorType* delta_x,
    std::size_t block_size)
{
    DenseVectorType const& F = vector_f;
    SparseMatrixType Jt = J.transpose();
    SparseMatrixType H = Jt.multiply(J);
    SparseMatrixType H_diag = H.diagonal_matrix();

    /* Compute RHS. */
    DenseVectorType g = Jt.multiply(F);
    g.negate_self();

    /* Add regularization to H. */
    H.mult_diagonal(1.0 + 1.0 / this->opts.trust_region_radius);

    Status status;
    if (block_size == 0)
    {
        /* Use preconditioned CG using the diagonal of H. */
        SparseMatrixType precond = H.diagonal_matrix();
        precond.cwise_invert();

        typedef sfm::ba::ConjugateGradient<double> CGSolver;
        CGSolver::Options cg_opts;
        cg_opts.max_iterations = this->opts.cg_max_iterations;
        cg_opts.tolerance = 1e-20;
        CGSolver solver(cg_opts);
        CGSolver::Status cg_status;
        cg_status = solver.solve(H, g, delta_x, &precond);
        status.num_cg_iterations = cg_status.num_iterations;

        switch (cg_status.info)
        {
            case CGSolver::CG_CONVERGENCE:
                status.success = true;
                break;
            case CGSolver::CG_MAX_ITERATIONS:
                status.success = true;
                break;
            case CGSolver::CG_INVALID_INPUT:
                std::cout << "BA: CG failed (invalid input)" << std::endl;
                status.success = false;
                return status;
            default:
                break;
        }
    }
    else if (block_size == 3)
    {
        /* Invert blocks of H directly */
        invert_block_matrix_3x3_inplace(&H);
        *delta_x = H.multiply(g);
        status.success = true;
        status.num_cg_iterations = 0;
    }
    else
    {
        status.success = false;
        throw std::invalid_argument("Unsupported block_size in linear solver");
    }

    status.predicted_error_decrease = delta_x->dot(H_diag.multiply(
        *delta_x).multiply(1.0 / this->opts.trust_region_radius).add(g));

    return status;
}

SFM_BA_NAMESPACE_END
SFM_NAMESPACE_END
