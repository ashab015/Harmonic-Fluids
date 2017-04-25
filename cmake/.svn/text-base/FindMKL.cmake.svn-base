###############################################################################
#  File: FindMKL.cmake
#  Copyright (c) 2008 by Changxi Zheng
# 
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

# - Find Intel MKL
# set variables
#
#  MKL_INCLUDE_DIR      MKL include directory
#  MKL_LIB_DIRS         MKL library directory
#
FIND_PATH(MKL_INCLUDE_DIR mkl.h
    PATHS ENV INCLUDE
    PATHS ${SYSTEM_INC_PATH})

# Find the mkl core library
FIND_LIBRARY(_mkl_lib_dir mkl_core
    PATHS ENV LD_LIBRARY_PATH
    PATHS ${SYSTEM_LIB_PATH}
    PATHS ENV LIBRARY_PATH)

# Find guide for MKL in intel icc 11 or later it may be
# located in different directory from where mkl_core is
FIND_LIBRARY(_mkl_guide_lib_dir guide
    PATHS ENV LD_LIBRARY_PATH
    PATHS ${SYSTEM_LIB_PATH}
    PATHS ENV LIBRARY_PATH)

IF (MKL_INCLUDE_DIR AND _mkl_lib_dir AND _mkl_guide_lib_dir)
    IF (NOT MKL_FIND_QUIETLY)
        MESSAGE(STATUS "Found MKL at: ${_mkl_lib_dir}")
    ENDIF (NOT MKL_FIND_QUIETLY)
    # encapsulate MKL_LIB_DIRS
    GET_FILENAME_COMPONENT(MKL_LIB_DIRS "${_mkl_lib_dir}" PATH)

    GET_FILENAME_COMPONENT(_mkl_lib_dir "${_mkl_guide_lib_dir}" PATH)
    LIST(APPEND MKL_LIB_DIRS ${_mkl_lib_dir})
    LIST(REMOVE_DUPLICATES MKL_LIB_DIRS)
ELSE (MKL_INCLUDE_DIR AND _mkl_lib_dir AND _mkl_guide_lib_dir)
    # Cannot find the library
    IF (MKL_FIND_REQUIRED)
        MESSAGE(FATAL_ERROR "Could not find MKL: $ENV{INCLUDE} $ENV{LD_LIBRARY_PATH}")
    ELSE (MKL_FIND_REQUIRED)
        IF ( NOT _mkl_lib_dir )
            MESSAGE(STATUS "WARNING: Could not find MKL library: $ENV{LD_LIBRARY_PATH}")
        ELSEIF (NOT _mkl_guide_lib_dir )
            MESSAGE(STATUS "WARNING: Could not find MKL guide library: $ENV{LD_LIBRARY_PATH}")
        ELSE ( NOT _mkl_lib_dir )
            MESSAGE(STATUS "WARNING: Could not find MKL include: $ENV{INCLUDE}")
        ENDIF ( NOT _mkl_lib_dir )
    ENDIF (MKL_FIND_REQUIRED)
ENDIF (MKL_INCLUDE_DIR AND _mkl_lib_dir AND _mkl_guide_lib_dir)
