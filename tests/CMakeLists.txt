set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")

include_directories("${PROJET_SOURCE_DIR}/src")

add_executable(runTests polyATest.cpp)
target_link_libraries(runTests aludetect)
