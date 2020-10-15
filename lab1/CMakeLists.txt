project(lab1)

add_executable(lab1 lab1.cpp)

target_link_libraries(lab1
    PRIVATE
        ${MPI_CXX_LIBRARIES})