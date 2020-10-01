project(lab1)

include_directories(/usr/include/mpi/)
set(CMAKE_CXX_COMPILER "/usr/bin/mpic++")

add_executable(lab1 lab1.cpp)