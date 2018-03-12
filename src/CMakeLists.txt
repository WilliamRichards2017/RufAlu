# =================================
# AluDetect
#
# src/CMakeLists.txt
# =================================

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/src/)

# set header and source files
set(ALUDETECT_SRC
	aluDetect.cpp
)

ADD_DEFINITIONS(-DBOOST_FALLTHROUGH)
INCLUDE_DIRECTORIES(
   ${BAMTOOLS_INCLUDE}
 )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") # clang Doesnt use pthread
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()

set(CMAKE_BUILD_TYPE Debug) # uncomment this when you need a debug build


#create binary
add_executable(aluDetect
  ${ALUDETECT_SRC}
)

# link libraries to executables
target_link_libraries(aluDetect
	${CORE_LIB}
)




