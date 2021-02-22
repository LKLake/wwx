cmake_minimum_required(VERSION 3.17)
project(wwx)

set(CMAKE_CXX_STANDARD 98)
include_directories("/usr/include/eigen3")

add_executable(wwx C-V.cpp)
target_link_libraries(wwx umfpack)