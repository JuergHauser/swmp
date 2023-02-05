
cmake_minimum_required(VERSION 3.11.0)
enable_language( Fortran )
project(swmp)
add_subdirectory(modules)
add_subdirectory(apps)

file (COPY ${CMAKE_CURRENT_SOURCE_DIR}/demos DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file (COPY ${CMAKE_CURRENT_SOURCE_DIR}/python DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

