/*
 * Copyright (C) 2015, Simon Fuhrmann
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include <iostream>
#include <fstream>
#include <stdexcept>

#include "util/timer.h"
#include "math/functions.h"
#include "math/matrix.h"
#include "math/matrix_tools.h"
#include "core/image_io.h"
#include "core/image_tools.h"
#include "features/sift.h"

FEATURES_NAMESPACE_BEGIN

Sift::Sift (Options const& options)
    : options(options)
{
    if (this->options.min_octave < -1
        || this->options.min_octave > this->options.max_octave)
        throw std::invalid_argument("Invalid octave range");

    if (this->options.contrast_threshold < 0.0f)
        this->options.contrast_threshold = 0.02f
            / static_cast<float>(this->options.num_samples_per_octave);

    if (this->options.debug_output)
        this->options.verbose_output = true;
}

/* ---------------------------------------------------------------- */

void
Sift::process (void)
{
    util::ClockTimer timer, total_timer;

    /*
     * Creates the scale space representation of the image by
     * sampling the scale space and computing the DoG images.
     * See Section 3, 3.2 and 3.3 in SIFT article.
     */
    if (this->options.verbose_output)
    {
        std::cout << "SIFT: Creating "
            << (this->options.max_octave - this->options.min_octave)
            << " octaves (" << this->options.min_octave << " to "
            << this->options.max_octave << ")..." << std::endl;
    }
    timer.reset();
    this->create_octaves();
    if (this->options.debug_output)
    {
        std::cout << "SIFT: Creating octaves took "
            << timer.get_elapsed() << "ms." << std::endl;
    }

    /*
     * Detects local extrema in the DoG function as described in Section 3.1.
     */
    if (this->options.debug_output)
    {
        std::cout << "SIFT: Detecting local extrema..." << std::endl;
    }
    timer.reset();
    this->extrema_detection();
    if (this->options.debug_output)
    {
        std::cout << "SIFT: Detected " << this->keypoints.size()
            << " keypoints, took " << timer.get_elapsed() << "ms." << std::endl;
    }

    /*
     * Accurate keypoint localization and filtering.
     * According to Section 4 in SIFT article.
     */
    if (this->options.debug_output)
    {
        std::cout << "SIFT: Localizing and filtering keypoints..." << std::endl;
    }
    timer.reset();
    this->keypoint_localization();
    if (this->options.debug_output)
    {
        std::cout << "SIFT: Retained " << this->keypoints.size() << " stable "
            << "keypoints, took " << timer.get_elapsed() << "ms." << std::endl;
    }

    /*
     * Difference of Gaussian images are not needed anymore.
     */
    for (std::size_t i = 0; i < this->octaves.size(); ++i)
        this->octaves[i].dog.clear();

    /*
     * Generate the list of keypoint descriptors.
     * See Section 5 and 6 in the SIFT article.
     * This list can in general be larger than the amount of keypoints,
     * since for each keypoint several descriptors may be created.
     */
    if (this->options.verbose_output)
    {
        std::cout << "SIFT: Generating keypoint descriptors..." << std::endl;
    }
    timer.reset();
    this->descriptor_generation();
    if (this->options.debug_output)
    {
        std::cout << "SIFT: Generated " << this->descriptors.size()
            << " descriptors, took " << timer.get_elapsed() << "ms."
            << std::endl;
    }
    if (this->options.verbose_output)
    {
        std::cout << "SIFT: Generated " << this->descriptors.size()
            << " descriptors from " << this->keypoints.size() << " keypoints,"
            << " took " << total_timer.get_elapsed() << "ms." << std::endl;
    }

    /* Free memory. */
    this->octaves.clear();
}

/* ---------------------------------------------------------------- */

void
Sift::set_image (core::ByteImage::ConstPtr img)
{
    if (img->channels() != 1 && img->channels() != 3)
        throw std::invalid_argument("Gray or color image expected");

    // 将图像转化成灰度图
    this->orig = core::image::byte_to_float_image(img);
    if (img->channels() == 3) {
        this->orig = core::image::desaturate<float>
            (this->orig, core::image::DESATURATE_AVERAGE);
    }
}

/* ---------------------------------------------------------------- */

void
Sift::set_float_image (core::FloatImage::ConstPtr img)
{
    if (img->channels() != 1 && img->channels() != 3)
        throw std::invalid_argument("Gray or color image expected");

    if (img->channels() == 3)
    {
        this->orig = core::image::desaturate<float>
            (img, core::image::DESATURATE_AVERAGE);
    }
    else
    {
        this->orig = img->duplicate();
    }
}

/* ---------------------------------------------------------------- */

void
Sift::create_octaves (void)
{
    this->octaves.clear();

    /*
     * Create octave -1. The original image is assumed to have blur
     * sigma = 0.5. The double size image therefore has sigma = 1.
     */
    if (this->options.min_octave < 0)
    {
        //std::cout << "Creating octave -1..." << std::endl;
        core::FloatImage::Ptr img
            = core::image::rescale_double_size_supersample<float>(this->orig);
        this->add_octave(img, this->options.inherent_blur_sigma * 2.0f,
            this->options.base_blur_sigma);
    }

    /*
     * Prepare image for the first positive octave by downsampling.
     * This code is executed only if min_octave > 0.
     */
    core::FloatImage::ConstPtr img = this->orig;
    for (int i = 0; i < this->options.min_octave; ++i)
        img = core::image::rescale_half_size_gaussian<float>(img);

    /*
     * Create new octave from 'img', then subsample octave image where
     * sigma is doubled to get a new base image for the next octave.
     */
    float img_sigma = this->options.inherent_blur_sigma;
    for (int i = std::max(0, this->options.min_octave);
        i <= this->options.max_octave; ++i)
    {
        //std::cout << "Creating octave " << i << "..." << std::endl;
        this->add_octave(img, img_sigma, this->options.base_blur_sigma);

        core::FloatImage::ConstPtr pre_base = octaves[octaves.size()-1].img[0];
        img = core::image::rescale_half_size_gaussian<float>(pre_base);

        img_sigma = this->options.base_blur_sigma;
    }
}

/* ---------------------------------------------------------------- */

void
Sift::add_octave (core::FloatImage::ConstPtr image,
        float has_sigma, float target_sigma)
{
    /*
     * First, bring the provided image to the target blur.
     * Since L * g(sigma1) * g(sigma2) = L * g(sqrt(sigma1^2 + sigma2^2)),
     * we need to blur with sigma = sqrt(target_sigma^2 - has_sigma^2).
     */
    float sigma = std::sqrt(MATH_POW2(target_sigma) - MATH_POW2(has_sigma));
    //std::cout << "Pre-blurring image to sigma " << target_sigma << " (has "
    //    << has_sigma << ", blur = " << sigma << ")..." << std::endl;
    core::FloatImage::Ptr base = (target_sigma > has_sigma
        ? core::image::blur_gaussian<float>(image, sigma)
        : image->duplicate());

    /* Create the new octave and add initial image. */
    this->octaves.push_back(Octave());
    Octave& oct = this->octaves.back();
    oct.img.push_back(base);

    /* 'k' is the constant factor between the scales in scale space. */
    float const k = std::pow(2.0f, 1.0f / this->options.num_samples_per_octave);
    sigma = target_sigma;

    /* Create other (s+2) samples of the octave to get a total of (s+3). */
    for (int i = 1; i < this->options.num_samples_per_octave + 3; ++i)
    {
        /* Calculate the blur sigma the image will get. */
        float sigmak = sigma * k;
        float blur_sigma = std::sqrt(MATH_POW2(sigmak) - MATH_POW2(sigma));

        /* Blur the image to create a new scale space sample. */
        //std::cout << "Blurring image to sigma " << sigmak << " (has " << sigma
        //    << ", blur = " << blur_sigma << ")..." << std::endl;
        core::FloatImage::Ptr img = core::image::blur_gaussian<float>
            (base, blur_sigma);
        oct.img.push_back(img);

        /* Create the Difference of Gaussian image (DoG). */
        //计算差分拉普拉斯 // todo revised by sway
        core::FloatImage::Ptr dog = core::image::subtract<float>(img, base);
        oct.dog.push_back(dog);

        /* Update previous image and sigma for next round. */
        base = img;
        sigma = sigmak;
    }
}

/* ---------------------------------------------------------------- */

void
Sift::extrema_detection (void)
{
    /* Delete previous keypoints. */
    this->keypoints.clear();

    /* Detect keypoints in each octave... */
    for (std::size_t i = 0; i < this->octaves.size(); ++i)
    {
        Octave const& oct(this->octaves[i]);
        /* In each octave, take three subsequent DoG images and detect. */
        for (int s = 0; s < (int)oct.dog.size() - 2; ++s)
        {
            core::FloatImage::ConstPtr samples[3] =
            { oct.dog[s + 0], oct.dog[s + 1], oct.dog[s + 2] };
            this->extrema_detection(samples, static_cast<int>(i)
                + this->options.min_octave, s);
        }
    }
}

/* ---------------------------------------------------------------- */

std::size_t
Sift::extrema_detection (core::FloatImage::ConstPtr s[3], int oi, int si)
{
    int const w = s[1]->width();
    int const h = s[1]->height();

    /* Offsets for the 9-neighborhood w.r.t. center pixel. */
    int noff[9] = { -1 - w, 0 - w, 1 - w, -1, 0, 1, -1 + w, 0 + w, 1 + w };

    /*
     * Iterate over all pixels in s[1], and check if pixel is maximum
     * (or minumum) in its 27-neighborhood.
     */
    int detected = 0;
    int off = w;
    for (int y = 1; y < h - 1; ++y, off += w)
        for (int x = 1; x < w - 1; ++x)
        {
            int idx = off + x;

            bool largest = true;
            bool smallest = true;
            float center_value = s[1]->at(idx);
            for (int l = 0; (largest || smallest) && l < 3; ++l)
                for (int i = 0; (largest || smallest) && i < 9; ++i)
                {
                    if (l == 1 && i == 4) // Skip center pixel
                        continue;
                    if (s[l]->at(idx + noff[i]) >= center_value)
                        largest = false;
                    if (s[l]->at(idx + noff[i]) <= center_value)
                        smallest = false;
                }

            /* Skip non-maximum values. */
            if (!smallest && !largest)
                continue;

            /* Yummy. Add detected scale space extremum. */
            Keypoint kp;
            kp.octave = oi;
            kp.x = static_cast<float>(x);
            kp.y = static_cast<float>(y);
            kp.sample = static_cast<float>(si);
            this->keypoints.push_back(kp);
            detected += 1;
        }

    return detected;
}

/* ---------------------------------------------------------------- */

void
Sift::keypoint_localization (void)
{
    /*
     * Iterate over all keypoints, accurately localize minima and maxima
     * in the DoG function by fitting a quadratic Taylor polynomial
     * around the keypoint.
     */

    int num_singular = 0;
    int num_keypoints = 0; // Write iterator
    for (std::size_t i = 0; i < this->keypoints.size(); ++i)
    {
        /* Copy keypoint. */
        Keypoint kp(this->keypoints[i]);

        /* Get corresponding octave and DoG images. */
        Octave const& oct(this->octaves[kp.octave - this->options.min_octave]);
        int sample = static_cast<int>(kp.sample);
        core::FloatImage::ConstPtr dogs[3] = { oct.dog[sample + 0], oct.dog[sample + 1], oct.dog[sample + 2] };

        /* Shorthand for image width and height. */
        int const w = dogs[0]->width();
        int const h = dogs[0]->height();
        /* The integer and floating point location of the keypoints. */
        int ix = static_cast<int>(kp.x);
        int iy = static_cast<int>(kp.y);
        int is = static_cast<int>(kp.sample);
        float delta_x, delta_y, delta_s;
        /* The first and second order derivatives. */
        float Dx, Dy, Ds;
        float Dxx, Dyy, Dss;
        float Dxy, Dxs, Dys;

        /*
         * Locate the keypoint using second order Taylor approximation.
         * The procedure might get iterated around a neighboring pixel if
         * the accurate keypoint is off by >0.6 from the center pixel.
         */
#       define AT(S,OFF) (dogs[S]->at(px + OFF))
        for (int j = 0; j < 5; ++j)
        {
            std::size_t px = iy * w + ix;

            /* Compute first and second derivatives. */
            Dx = (AT(1,1) - AT(1,-1)) * 0.5f;
            Dy = (AT(1,w) - AT(1,-w)) * 0.5f;
            Ds = (AT(2,0) - AT(0,0))  * 0.5f;

            Dxx = AT(1,1) + AT(1,-1) - 2.0f * AT(1,0);
            Dyy = AT(1,w) + AT(1,-w) - 2.0f * AT(1,0);
            Dss = AT(2,0) + AT(0,0)  - 2.0f * AT(1,0);

            Dxy = (AT(1,1+w) + AT(1,-1-w) - AT(1,-1+w) - AT(1,1-w)) * 0.25f;
            Dxs = (AT(2,1)   + AT(0,-1)   - AT(2,-1)   - AT(0,1))   * 0.25f;
            Dys = (AT(2,w)   + AT(0,-w)   - AT(2,-w)   - AT(0,w))   * 0.25f;

            /* Setup the Hessian matrix. */
            math::Matrix3f H;
            /****************************task-1-0  构造Hessian矩阵 ******************************/
            /*
             * 参考第32页slide的Hessian矩阵构造方式填充H矩阵，其中dx=dy=d_sigma=1, 其中A矩阵按照行顺序存储，即
             * H=[H[0], H[1], H[2]]
             *   [H[3], H[4], H[5]]
             *   [H[6], H[7], H[8]]
             */

            /**********************************************************************************/
            H[0] = Dxx; H[1] = Dxy; H[2] = Dxs;
            H[3] = Dxy; H[4] = Dyy; H[5] = Dys;
            H[6] = Dxs; H[7] = Dys; H[8] = Dss;


            /* Compute determinant to detect singular matrix. */
            float detH = math::matrix_determinant(H);
            if (MATH_EPSILON_EQ(detH, 0.0f, 1e-15f))
            {
                num_singular += 1;
                delta_x = delta_y = delta_s = 0.0f; // FIXME: Handle this case?
                break;
            }
            /* Invert the matrix to get the accurate keypoint. */
            math::Matrix3f H_inv = math::matrix_inverse(H, detH);
            math::Vec3f b(-Dx, -Dy, -Ds);


            //math::Vec3f delta;
            /****************************task-1-1  求解偏移量deta ******************************/

             /* 参考第30页slide delta_x的求解方式 delta_x = inv(H)*b
             * 请在此处给出delta的表达式
             */
                     /*                  */
                     /*    此处添加代码    */
                     /*                  */

            /**********************************************************************************/


            math::Vec3f delta = H_inv * b;
            delta_x = delta[0];
            delta_y = delta[1];
            delta_s = delta[2];


            /* Check if accurate location is far away from pixel center. */
            // dx =0 表示|dx|>0.6f
            int dx = (delta_x > 0.6f && ix < w-2) * 1 + (delta_x < -0.6f && ix > 1) * -1;
            int dy = (delta_y > 0.6f && iy < h-2) * 1 + (delta_y < -0.6f && iy > 1) * -1;

            /* If the accurate location is closer to another pixel,
             * repeat localization around the other pixel. */
            if (dx != 0 || dy != 0)
            {
                ix += dx;
                iy += dy;
                continue;
            }
            /* Accurate location looks good. */
            break;
        }


        /* Calcualte function value D(x) at accurate keypoint x. */
        /*****************************task1-2求解极值点处的DoG值val ***************************/
         /*
          * 参考第30页slides的机极值点f(x)的求解公式f(x) = f(x0) + 0.5* delta.dot(D)
          * 其中
          * f(x0)--表示插值点(ix, iy, is) 处的DoG值，可通过dogs[1]->at(ix, iy, 0)获取
          * delta--为上述求得的delta=[delta_x, delta_y, delta_s]
          * D--为一阶导数，表示为(Dx, Dy, Ds)
          * 请给出求解val的代码
          */
        //float val = 0.0;
        /*                  */
        /*    此处添加代码    */
        /*                  */
        /************************************************************************************/
        float val = dogs[1]->at(ix, iy, 0) + 0.5f * (Dx * delta_x + Dy * delta_y + Ds * delta_s);
        /* Calcualte edge response score Tr(H)^2 / Det(H), see Section 4.1. */

         /**************************去除边缘点，参考第33页slide 仔细阅读代码 ****************************/
        float hessian_trace = Dxx + Dyy;
        float hessian_det = Dxx * Dyy - MATH_POW2(Dxy);
        float hessian_score = MATH_POW2(hessian_trace) / hessian_det;
        float score_thres = MATH_POW2(this->options.edge_ratio_threshold + 1.0f)
            / this->options.edge_ratio_threshold;
        /********************************************************************************/

        /*
         * Set accurate final keypoint location.
         */
        kp.x = (float)ix + delta_x;
        kp.y = (float)iy + delta_y;
        kp.sample = (float)is + delta_s;

        /*
         * Discard keypoints with:
         * 1. low contrast (value of DoG function at keypoint),
         * 2. negative hessian determinant (curvatures with different sign),
         *    Note that negative score implies negative determinant.
         * 3. large edge response (large hessian score),
         * 4. unstable keypoint accurate locations,
         * 5. keypoints beyond the scale space boundary.
         */
        if (std::abs(val) < this->options.contrast_threshold
            || hessian_score < 0.0f || hessian_score > score_thres
            || std::abs(delta_x) > 1.5f || std::abs(delta_y) > 1.5f || std::abs(delta_s) > 1.0f
            || kp.sample < -1.0f
            || kp.sample > (float)this->options.num_samples_per_octave
            || kp.x < 0.0f || kp.x > (float)(w - 1)
            || kp.y < 0.0f || kp.y > (float)(h - 1))
        {
            //std::cout << " REJECTED!" << std::endl;
            continue;
        }

        /* Keypoint is accepted, copy to write iter and advance. */
        this->keypoints[num_keypoints] = kp;
        num_keypoints += 1;
    }

    /* Limit vector size to number of accepted keypoints. */
    this->keypoints.resize(num_keypoints);

    if (this->options.debug_output && num_singular > 0)
    {
        std::cout << "SIFT: Warning: " << num_singular
            << " singular matrices detected!" << std::endl;
    }
}

/* ---------------------------------------------------------------- */

void
Sift::descriptor_generation (void)
{
    if (this->octaves.empty())
        throw std::runtime_error("Octaves not available!");
    if (this->keypoints.empty())
        return;

    this->descriptors.clear();
    this->descriptors.reserve(this->keypoints.size() * 3 / 2);

    /*
     * Keep a buffer of S+3 gradient and orientation images for the current
     * octave. Once the octave is changed, these images are recomputed.
     * To ensure efficiency, the octave index must always increase, never
     * decrease, which is enforced during the algorithm.
     */
    int octave_index = this->keypoints[0].octave;
    Octave* octave = &this->octaves[octave_index - this->options.min_octave];

    // todo 计算每个octave中所有图像的梯度值和方向，具体得, octave::grad存储图像的梯度响应值，octave::ori存储梯度方向
    this->generate_grad_ori_images(octave);

    /* Walk over all keypoints and compute descriptors. */
    for (std::size_t i = 0; i < this->keypoints.size(); ++i)
    {
        Keypoint const& kp(this->keypoints[i]);

        /* Generate new gradient and orientation images if octave changed. */
        if (kp.octave > octave_index)
        {
            /* Clear old octave gradient and orientation images. */
            if (octave)
            {
                octave->grad.clear();
                octave->ori.clear();
            }
            /* Setup new octave gradient and orientation images. */
            octave_index = kp.octave;
            octave = &this->octaves[octave_index - this->options.min_octave];
            this->generate_grad_ori_images(octave);
        }
        else if (kp.octave < octave_index)
        {
            throw std::runtime_error("Decreasing octave index!");
        }

        /* Orientation assignment. This returns multiple orientations. */
        /* todo 统计直方图找到特征点主方向,找到几个主方向*/
        std::vector<float> orientations;
        orientations.reserve(8);
        this->orientation_assignment(kp, octave, orientations);

        /* todo 生成特征向量,同一个特征点可能有多个描述子，为了提升匹配的稳定性*/
        /* Feature vector extraction. */
        for (std::size_t j = 0; j < orientations.size(); ++j)
        {
            Descriptor desc;
            float const scale_factor = std::pow(2.0f, kp.octave);
            desc.x = scale_factor * (kp.x + 0.5f) - 0.5f;
            desc.y = scale_factor * (kp.y + 0.5f) - 0.5f;
            desc.scale = this->keypoint_absolute_scale(kp);
            desc.orientation = orientations[j];
            if (this->descriptor_assignment(kp, desc, octave))
                this->descriptors.push_back(desc);
        }
    }
}

/* ---------------------------------------------------------------- */

void
Sift::generate_grad_ori_images (Octave* octave)
{
    octave->grad.clear();
    octave->grad.reserve(octave->img.size());
    octave->ori.clear();
    octave->ori.reserve(octave->img.size());

    int const width = octave->img[0]->width();
    int const height = octave->img[0]->height();

    //std::cout << "Generating gradient and orientation images..." << std::endl;
    for (std::size_t i = 0; i < octave->img.size(); ++i)
    {
        core::FloatImage::ConstPtr img = octave->img[i];
        core::FloatImage::Ptr grad = core::FloatImage::create(width, height, 1);
        core::FloatImage::Ptr ori = core::FloatImage::create(width, height, 1);

        int image_iter = width + 1;
        for (int y = 1; y < height - 1; ++y, image_iter += 2)
            for (int x = 1; x < width - 1; ++x, ++image_iter)
            {
                float m1x = img->at(image_iter - 1);
                float p1x = img->at(image_iter + 1);
                float m1y = img->at(image_iter - width);
                float p1y = img->at(image_iter + width);
                float dx = 0.5f * (p1x - m1x);
                float dy = 0.5f * (p1y - m1y);

                float atan2f = std::atan2(dy, dx);
                grad->at(image_iter) = std::sqrt(dx * dx + dy * dy);
                ori->at(image_iter) = atan2f < 0.0f
                    ? atan2f + MATH_PI * 2.0f : atan2f;
            }
        octave->grad.push_back(grad);
        octave->ori.push_back(ori);
    }
}

/* ---------------------------------------------------------------- */

void
Sift::orientation_assignment (Keypoint const& kp,
    Octave const* octave, std::vector<float>& orientations)
{
    int const nbins = 36;
    float const nbinsf = static_cast<float>(nbins);

    /* Prepare 36-bin histogram. */
    float hist[nbins];
    std::fill(hist, hist + nbins, 0.0f);

    /* Integral x and y coordinates and closest scale sample. */
    int const ix = static_cast<int>(kp.x + 0.5f);
    int const iy = static_cast<int>(kp.y + 0.5f);
    int const is = static_cast<int>(math::round(kp.sample));
    float const sigma = this->keypoint_relative_scale(kp);

    /* Images with its dimension for the keypoint. */
    core::FloatImage::ConstPtr grad(octave->grad[is + 1]);
    core::FloatImage::ConstPtr ori(octave->ori[is + 1]);
    int const width = grad->width();
    int const height = grad->height();

    /*
     * Compute window size 'win', the full window has  2 * win + 1  pixel.
     * The factor 3 makes the window large enough such that the gaussian
     * has very little weight beyond the window. The value 1.5 is from
     * the SIFT paper. If the window goes beyond the image boundaries,
     * the keypoint is discarded.
     */
    float const sigma_factor = 1.5f;
    int win = static_cast<int>(sigma * sigma_factor * 3.0f);
    if (ix < win || ix + win >= width || iy < win || iy + win >= height)
        return;

    /* Center of keypoint index. */
    int center = iy * width + ix;
    float const dxf = kp.x - static_cast<float>(ix);
    float const dyf = kp.y - static_cast<float>(iy);
    float const maxdist = static_cast<float>(win*win) + 0.5f;

    /* Populate histogram over window, intersected with (1,1), (w-2,h-2). */
    for (int dy = -win; dy <= win; ++dy)
    {
        int const yoff = dy * width;
        for (int dx = -win; dx <= win; ++dx)
        {
            /* Limit to circular window (centered at accurate keypoint). */
            float const dist = MATH_POW2(dx-dxf) + MATH_POW2(dy-dyf);
            if (dist > maxdist)
                continue;

            float gm = grad->at(center + yoff + dx); // gradient magnitude
            float go = ori->at(center + yoff + dx); // gradient orientation
            float weight = math::gaussian_xx(dist, sigma * sigma_factor);
            int bin = static_cast<int>(nbinsf * go / (2.0f * MATH_PI));
            bin = math::clamp(bin, 0, nbins - 1);
            hist[bin] += gm * weight;
        }
    }

    /* Smooth histogram. */
    for (int i = 0; i < 6; ++i)
    {
        float first = hist[0];
        float prev = hist[nbins - 1];
        for (int j = 0; j < nbins - 1; ++j)
        {
            float current = hist[j];
            hist[j] = (prev + current + hist[j + 1]) / 3.0f;
            prev = current;
        }
        hist[nbins - 1] = (prev + hist[nbins - 1] + first) / 3.0f;
    }

    /* Find maximum element. */
    float maxh = *std::max_element(hist, hist + nbins);

    /* Find peaks within 80% of max element. */
    for (int i = 0; i < nbins; ++i)
    {
        float h0 = hist[(i + nbins - 1) % nbins];
        float h1 = hist[i];
        float h2 = hist[(i + 1) % nbins];

        /* These peaks must be a local maximum! */
        if (h1 <= 0.8f * maxh || h1 <= h0 || h1 <= h2)
            continue;

        /*
         * Quadratic interpolation to find accurate maximum.
         * f(x) = ax^2 + bx + c, f(-1) = h0, f(0) = h1, f(1) = h2
         * --> a = 1/2 (h0 - 2h1 + h2), b = 1/2 (h2 - h0), c = h1.
         * x = f'(x) = 2ax + b = 0 --> x = -1/2 * (h2 - h0) / (h0 - 2h1 + h2)
         */
        float x = -0.5f * (h2 - h0) / (h0 - 2.0f * h1 + h2);
        float o =  2.0f * MATH_PI * (x + (float)i + 0.5f) / nbinsf;
        orientations.push_back(o);
    }
}

/* ---------------------------------------------------------------- */

bool
Sift::descriptor_assignment (Keypoint const& kp, Descriptor& desc,
    Octave const* octave)
{
    /*
     * The final feature vector has size PXB * PXB * OHB.
     * The following constants should not be changed yet, as the
     * (PXB^2 * OHB = 128) element feature vector is still hard-coded.
     */
    //int const PIX = 16; // Descriptor region with 16x16 pixel
    int const PXB = 4; // Pixel bins with 4x4 bins
    int const OHB = 8; // Orientation histogram with 8 bins

    /* Integral x and y coordinates and closest scale sample. */
    int const ix = static_cast<int>(kp.x + 0.5f);
    int const iy = static_cast<int>(kp.y + 0.5f);
    int const is = static_cast<int>(math::round(kp.sample));
    float const dxf = kp.x - static_cast<float>(ix);
    float const dyf = kp.y - static_cast<float>(iy);
    float const sigma = this->keypoint_relative_scale(kp);

    /* Images with its dimension for the keypoint. */
    core::FloatImage::ConstPtr grad(octave->grad[is + 1]);
    core::FloatImage::ConstPtr ori(octave->ori[is + 1]);
    int const width = grad->width();
    int const height = grad->height();

    /* Clear feature vector. */
    desc.data.fill(0.0f);

    /* Rotation constants given by descriptor orientation. */
    float const sino = std::sin(desc.orientation);
    float const coso = std::cos(desc.orientation);

    /*
     * Compute window size.
     * Each spacial bin has an extension of 3 * sigma (sigma is the scale
     * of the keypoint). For interpolation we need another half bin at
     * both ends in each dimension. And since the window can be arbitrarily
     * rotated, we need to multiply with sqrt(2). The window size is:
     * 2W = sqrt(2) * 3 * sigma * (PXB + 1).
     */
    float const binsize = 3.0f * sigma;
    int win = MATH_SQRT2 * binsize * (float)(PXB + 1) * 0.5f;
    if (ix < win || ix + win >= width || iy < win || iy + win >= height)
        return false;

    /*
     * Iterate over the window, intersected with the image region
     * from (1,1) to (w-2, h-2) since gradients/orientations are
     * not defined at the boundary pixels. Add all samples to the
     * corresponding bin.
     */
    int const center = iy * width + ix; // Center pixel at KP location
    for (int dy = -win; dy <= win; ++dy)
    {
        int const yoff = dy * width;
        for (int dx = -win; dx <= win; ++dx)
        {
            /* Get pixel gradient magnitude and orientation. */
            float const mod = grad->at(center + yoff + dx);
            float const angle = ori->at(center + yoff + dx);
            float theta = angle - desc.orientation;
            if (theta < 0.0f)
                theta += 2.0f * MATH_PI;

            /* Compute fractional coordinates w.r.t. the window. */
            float const winx = (float)dx - dxf;
            float const winy = (float)dy - dyf;

            /*
             * Compute normalized coordinates w.r.t. bins. The window
             * coordinates are rotated around the keypoint. The bins are
             * chosen such that 0 is the coordinate of the first bins center
             * in each dimension. In other words, (0,0,0) is the coordinate
             * of the first bin center in the three dimensional histogram.
             */
            float binoff = (float)(PXB - 1) / 2.0f;
            float binx = (coso * winx + sino * winy) / binsize + binoff;
            float biny = (-sino * winx + coso * winy) / binsize + binoff;
            float bint = theta * (float)OHB / (2.0f * MATH_PI) - 0.5f;

            /* Compute circular window weight for the sample. */
            float gaussian_sigma = 0.5f * (float)PXB;
            float gaussian_weight = math::gaussian_xx
                (MATH_POW2(binx - binoff) + MATH_POW2(biny - binoff),
                gaussian_sigma);

            /* Total contribution of the sample in the histogram is now: */
            float contrib = mod * gaussian_weight;

            /*
             * Distribute values into bins (using trilinear interpolation).
             * Each sample is inserted into 8 bins. Some of these bins may
             * not exist, because the sample is outside the keypoint window.
             */
            int bxi[2] = { (int)std::floor(binx), (int)std::floor(binx) + 1 };
            int byi[2] = { (int)std::floor(biny), (int)std::floor(biny) + 1 };
            int bti[2] = { (int)std::floor(bint), (int)std::floor(bint) + 1 };

            float weights[3][2] = {
                { (float)bxi[1] - binx, 1.0f - ((float)bxi[1] - binx) },
                { (float)byi[1] - biny, 1.0f - ((float)byi[1] - biny) },
                { (float)bti[1] - bint, 1.0f - ((float)bti[1] - bint) }
            };

            // Wrap around orientation histogram
            if (bti[0] < 0)
                bti[0] += OHB;
            if (bti[1] >= OHB)
                bti[1] -= OHB;

            /* Iterate the 8 bins and add weighted contrib to each. */
            int const xstride = OHB;
            int const ystride = OHB * PXB;
            for (int y = 0; y < 2; ++y)
                for (int x = 0; x < 2; ++x)
                    for (int t = 0; t < 2; ++t)
                    {
                        if (bxi[x] < 0 || bxi[x] >= PXB
                            || byi[y] < 0 || byi[y] >= PXB)
                            continue;

                        int idx = bti[t] + bxi[x] * xstride + byi[y] * ystride;
                        desc.data[idx] += contrib * weights[0][x]
                            * weights[1][y] * weights[2][t];
                    }
        }
    }

    /* Normalize the feature vector. */
    desc.data.normalize();

    /* Truncate descriptor values to 0.2. */
    for (int i = 0; i < PXB * PXB * OHB; ++i)
        desc.data[i] = std::min(desc.data[i], 0.2f);

    /* Normalize once again. */
    desc.data.normalize();

    return true;
}

/* ---------------------------------------------------------------- */

/*
 * The scale of a keypoint is: scale = sigma0 * 2^(octave + (s+1)/S).
 * sigma0 is the initial blur (1.6), octave the octave index of the
 * keypoint (-1, 0, 1, ...) and scale space sample s in [-1,S+1] where
 * S is the amount of samples per octave. Since the initial blur 1.6
 * corresponds to scale space sample -1, we add 1 to the scale index.
 */

float
Sift::keypoint_relative_scale (Keypoint const& kp)
{
    return this->options.base_blur_sigma * std::pow(2.0f,
        (kp.sample + 1.0f) / this->options.num_samples_per_octave);
}

float
Sift::keypoint_absolute_scale (Keypoint const& kp)
{
    return this->options.base_blur_sigma * std::pow(2.0f,
        kp.octave + (kp.sample + 1.0f) / this->options.num_samples_per_octave);
}

/* ---------------------------------------------------------------- */

void
Sift::load_lowe_descriptors (std::string const& filename, Descriptors* result)
{
    std::ifstream in(filename.c_str());
    if (!in.good())
        throw std::runtime_error("Cannot open descriptor file");

    int num_descriptors;
    int num_dimensions;
    in >> num_descriptors >> num_dimensions;
    if (num_descriptors > 100000 || num_dimensions != 128)
    {
        in.close();
        throw std::runtime_error("Invalid number of descriptors/dimensions");
    }
    result->clear();
    result->reserve(num_descriptors);
    for (int i = 0; i < num_descriptors; ++i)
    {
        Sift::Descriptor descriptor;
        in >> descriptor.y >> descriptor.x
            >> descriptor.scale >> descriptor.orientation;
        for (int j = 0; j < 128; ++j)
            in >> descriptor.data[j];
        descriptor.data.normalize();
        result->push_back(descriptor);
    }

    if (!in.good())
    {
        result->clear();
        in.close();
        throw std::runtime_error("Error while reading descriptors");
    }

    in.close();
}

FEATURES_NAMESPACE_END
