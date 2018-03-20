#  For more information, please see: http://software.sci.utah.edu
# 
#  The MIT License
# 
#  Copyright (c) 2015 Scientific Computing and Imaging Institute,
#  University of Utah.
# 
#  
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
# 
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software. 
# 
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.

SET_PROPERTY(DIRECTORY PROPERTY "EP_BASE" ${ep_base})


# Setting up external library for Minimap2
SET(MINIMAP2_PROJECT minimap2_project CACHE INTERNAL "minimap2 project name")
SET(MINIMAP2_DIR ${CMAKE_BINARY_DIR}/externals/minimap2 CACHE INTERNAL "minimap2 project directory")
ExternalProject_Add(${MINIMAP2_PROJECT}
	GIT_REPOSITORY https://github.com/lh3/minimap2.git
	GIT_TAG master
        CONFIGURE_COMMAND ""
	BUILD_COMMAND "make"
	INSTALL_COMMAND ""
	UPDATE_COMMAND ""
	BUILD_IN_SOURCE 1
    PREFIX ${MINIMAP2_DIR}
    CMAKE_CACHE_ARGS
        -DCMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER:STRING=${CMAKE_CXX_COMPILER}

)

ExternalProject_Get_Property(${MINIMAP2_PROJECT} SOURCE_DIR)
ExternalProject_Get_Property(${MINIMAP2_PROJECT} BINARY_DIR)

MESSAGE("BINARY_DIR: ${BINARY_DIR}")

MESSAGE("SRC_DIR: ${SOURCE_DIR}")


SET(MINIMAP2_LIB ${BINARY_DIR}/libminimap2.a CACHE INTERNAL "MINIMAP2 Library")
SET(MINIMAP2_INCLUDE ${SOURCE_DIR} CACHE INTERNAL "MINIMAP2 Include")