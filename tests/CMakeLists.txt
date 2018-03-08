set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/)

INCLUDE_DIRECTORIES(
   ${BAMTOOLS_INCLUDE}
 )

add_executable(runTests tests.cpp)
target_link_libraries(runTests ${CORE_LIB})
