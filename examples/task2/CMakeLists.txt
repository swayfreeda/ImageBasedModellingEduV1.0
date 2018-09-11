project(class2)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-fPIC")

include_directories("../..")


#test math basic
set(MATH_BASIC_FILE
        class2_test_math_basic.cc
        )
add_executable(task2_test_math_basic ${MATH_BASIC_FILE} )
target_link_libraries(task2_test_math_basic  util core )


# test camera model
set(CAMERA_MODEL
        class2_test_camera.cc
        )

add_executable(task2_test_camera_model ${CAMERA_MODEL}  )
target_link_libraries(task2_test_camera_model util core features )


# test fundamental 8 point
set(FUNDAMENTAL_8_POINT_FILE
        class2_test_fundamental_8_point.cc)
add_executable(task2_test_fundamental_8_point ${FUNDAMENTAL_8_POINT_FILE} )
target_link_libraries(task2_test_fundamental_8_point sfm util core features )



# test fundamental ransac
set(FUNDAMENTAL_RANSAC_FILE
        class2_test_fundamental_ransac.cc)
add_executable(task2_test_fundamental_ransac ${FUNDAMENTAL_RANSAC_FILE})
target_link_libraries(task2_test_fundamental_ransac sfm util core features )


set(POSE_FROM_FUNDAMENTAL
        calss2_test_pose_from_fundamental.cc
        )
add_executable(task2_test_pose_from_fundamental ${POSE_FROM_FUNDAMENTAL} )
target_link_libraries(task2_test_pose_from_fundamental sfm util core features )
