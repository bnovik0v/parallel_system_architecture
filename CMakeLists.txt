cmake_minimum_required(VERSION 3.17)
project(parallel_system_architecture)

set(CMAKE_CXX_STANDARD 14)

add_executable(parallel_system_architecture main.cpp)

include_directories(lab1)
add_subdirectory(lab1)