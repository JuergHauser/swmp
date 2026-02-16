cmake_minimum_required(VERSION 3.15...3.26)

project(pyswmp LANGUAGES Fortran)

# Set verbose makefile for debugging
set(CMAKE_VERBOSE_MAKEFILE on)

# Compiler-specific flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    # GNU Fortran (gfortran) flags
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -fbounds-check -fbacktrace -Wall -Wextra -fdiagnostics-color=auto")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -ffast-math")
    set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O2 -g -fbacktrace")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    # Intel Fortran flags
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -check bounds -traceback -warn all")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -xHost -ipo")
    set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O2 -g -traceback")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI|NVHPC")
    # NVIDIA HPC SDK (formerly PGI) flags
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -C -traceback")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fast")
    set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O2 -g -traceback")

else()
    # Generic flags for unknown compilers
    message(WARNING "Unknown Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}. Using generic flags.")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
endif()

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type (Debug, Release, RelWithDebInfo)" FORCE)
    message(STATUS "Build type not specified, defaulting to Release")
endif()

message(STATUS "Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Fortran flags: ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}}")

add_subdirectory(swmp)

# Install C header
install(FILES swmp/swmp.h DESTINATION include/swmp)
