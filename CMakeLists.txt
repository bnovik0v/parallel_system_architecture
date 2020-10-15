cmake_minimum_required(VERSION 3.16)
project(parallel_system_architecture)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

message(STATUS "Configuring MPI")
find_package(MPI)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

message(STATUS "Configuring OpenMP")
find_package(OpenMP)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

add_subdirectory(lab1)
add_subdirectory(lab2)
add_subdirectory(lab3)
add_subdirectory(lab4)
add_subdirectory(lab5)