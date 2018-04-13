set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/)

INCLUDE_DIRECTORIES(
   ${BAMTOOLS_INCLUDE}
   ${ZLIB_INCLUDE}
   ${MINIMAP2_INCLUDE}
   ${GTEST_SOURCE_DIR}/include
   ${GTEST_SOURCE_DIR}
 )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") # clang Doesnt use pthread
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()

add_executable(rufalu_tests 
	gtest.cpp
)

target_link_libraries(rufalu_tests 
	${CORE_LIB}
	${GTEST_LIB}
)

add_dependencies(rufalu_tests ${ALUDETECT_EXTERNAL_PROJECT})
