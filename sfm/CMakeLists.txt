project(sfm)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-fPIC")


include_directories("..")
set(HEADERS
       defines.h
        camera_pose.h
        camera_database.h
        bundler_common.h
        feature_set.h
        ransac.h
        fundamental.h
        ransac_homography.h
        homography.h
        ransac_homography.h
        pose_p3p.h
        ransac_pose_p3p.h
        bundle_adjustment.h
        ba_types.h
        ba_linear_solver.h
        ba_sparse_matrix.h
        ba_dense_vector.h
        ba_conjugate_gradient.h
        ba_cholesky.h
        extract_focal_length.h
        triangulate.h
        )

set(SOURCE_FILES
        camera_database.cc
        bundler_common.cc
        feature_set.cc
        ransac.cc
        fundamental.cc
        ransac_homography.cc
        homography.cc
        ransac_fundamental.cc
        pose_p3p.cc
        ransac_pose_p3p.cc
        bundle_adjustment.cc
        ba_linear_solver.cc
        extract_focal_length.cc
        triangulate.cc
        )
add_library(sfm ${HEADERS} ${SOURCE_FILES})
target_link_libraries(sfm core util features)

