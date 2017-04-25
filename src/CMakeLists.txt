## CMake specification for src directory

set(SOURCES
    TerminalColor.h     Bubble.hpp
    TinyBubble.hpp      Constants.hpp
    HelmholtzSolutions.hpp 
    PrecisionDef.h
    SpecialFunc.hpp)

# cluster/
set(CLUSTER_SRC
    cluster/harmonic_fluids.hpp
    cluster/harmonic_fluids_svc.cc
    cluster/harmonic_fluids_xdr.cc
    cluster/MultipoleSolver.hpp
    cluster/MultipoleSolver.cpp)
# 

add_definitions(-Wall)
if (CMAKE_CXX_COMPILER MATCHES ".*icpc$")
    # if use intel's compiler, disable some warnings
    add_definitions(-wd981 -wd383 -wd444 -fast)
else (CMAKE_CXX_COMPILER MATCHES ".*icpc$")
    if ( USE_OPENMP )
        set(OPENMP_LIBS gomp)
    endif (USE_OPENMP )
endif (CMAKE_CXX_COMPILER MATCHES ".*icpc$")

include_directories(
        ${CMAKE_CURRENT_SOURCE_DIR} 
        ${CMAKE_CURRENT_BINARY_DIR}
        ${Boost_INCLUDE_DIR})
list(APPEND MATH_LIBS ${MKL_LIBS})
link_directories(${MKL_LIB_DIRS})

add_executable(hfsvc
    ${SOURCES}      ${CLUSTER_SRC})
target_link_libraries(hfsvc
    ${MKL_LIBS}     ${OPENMP_LIBS}
    ${CMAKE_THREAD_LIBS_INIT})
