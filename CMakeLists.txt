project(CircleFit)

cmake_minimum_required(VERSION 2.8.9)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_executable(CircleFit CircleFit.cxx )
target_link_libraries(CircleFit  ${ITK_LIBRARIES})
