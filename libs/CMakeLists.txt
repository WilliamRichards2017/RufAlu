MESSAGE("BAMTOOLS LIB DIR: ${BAMTOOLS_INCLUDE}")

INCLUDE_DIRECTORIES(
    ${BAMTOOLS_INCLUDE}
    ${ZLIB_INCLUDE}
    ${MINIMAP2_INCLUDE}
)

set(CMAKE_BUILD_TYPE Debug) # uncomment this when you need a debug build

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") # clang Doesnt use pthread
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")

add_library(aludetect_core STATIC
	polyATail.cpp
	discordantAndChimericReads.cpp
	knownAlus.cpp
	
  )

SET(CORE_LIB aludetect_core CACHE STRING "THE CORE LIBRARY")

TARGET_LINK_LIBRARIES(${CORE_LIB}
  ${BAMTOOLS_LIB}
  ${ZLIB_LIBRARY}
  ${MINIMAP2_LIB}
)

add_dependencies(${CORE_LIB} ${ALUDETECT_EXTERNAL_PROJECT})
