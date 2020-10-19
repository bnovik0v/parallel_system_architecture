project(lab3)

add_executable(lab3 lab3.cpp)

target_link_libraries(lab3
    PRIVATE
        ${MPI_CXX_LIBRARIES})