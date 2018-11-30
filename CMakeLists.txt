# =================================
# RufAlu
#
# RufALu/CMakeLists.txt
# =================================

#Set up min CMake version
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(RufAlu)

SET(VERSION_MAJOR "0.8")
SET(VERSION_MINOR "0")
SET(VERSION_PATCH "0")
SET(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

INCLUDE_DIRECTORIES(src, externals, libs)

UNSET(CORE_LIBS CACHE)
SET(CORE_LIBS CACHE LIST "List of plugin libs")

ADD_SUBDIRECTORY(externals)
ADD_SUBDIRECTORY(libs)
ADD_SUBDIRECTORY(src)


