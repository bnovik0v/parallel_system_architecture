project(lab2)

add_executable(lab2 lab2.cpp)

target_link_libraries(lab2
    PUBLIC
        OpenMP::OpenMP_CXX)