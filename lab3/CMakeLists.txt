project(lab3)

include_directories(/usr/include/mpi/)
set(CMAKE_CXX_COMPILER "/usr/bin/mpic++")

add_executable(lab3 lab3.cpp)