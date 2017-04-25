# CMake Project file for bubble radiation solver for Harmonic Fluids project
# Changxi Zheng (cxzheng@cs.cornell.edu)
# Copyright @ Jan 2009
project(hfsvc)

cmake_minimum_required(VERSION 2.6)

# General Configuration: find all the required libraries.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# compile option
option(USE_DEBUG "Turn on the debug mode" OFF)
option(USE_64BIT_ARCH "Compile the 64bit executable" ON)
option(USE_OPENMP "Turn on the OpenMP feature when compiling" ON)

#===================================================================
## Compiler
# set compiler flags for debug/release
if ( USE_DEBUG )
    add_definitions(-DDEBUG)
    set(CMAKE_BUILD_TYPE Debug)
else ( USE_DEBUG )
    set(CMAKE_BUILD_TYPE Release)
endif ( USE_DEBUG )

if ( USE_64BIT_ARCH )
    add_definitions(-m64)
    set(LINK_FLAGS -m64)
endif ( USE_64BIT_ARCH )

if ( USE_OPENMP )
    add_definitions(-DUSE_OPENMP)
endif ( USE_OPENMP )

#===================================================================
# Libraries
# check Intel MKL library
find_package(MKL REQUIRED)
if (MKL_INCLUDE_DIR AND MKL_LIB_DIRS) 
    set(MKL_LIBS guide mkl mkl_lapack mkl_core)
    # check the pthread
    find_package(Threads 1.0 REQUIRED)

    add_definitions(-DUSE_MKL)
endif ()

find_package(Boost 1.37 REQUIRED)
if ( Boost_FOUND )
    add_definitions(-DHAVE_BOOST)
endif () 

if ( USE_OPENMP )
    find_package(OpenMP REQUIRED)
    add_definitions(${OpenMP_CXX_FLAGS})
endif (USE_OPENMP )

add_subdirectory(src)

