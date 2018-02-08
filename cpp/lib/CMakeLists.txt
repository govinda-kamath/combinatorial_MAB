cmake_minimum_required(VERSION 3.0)

include_directories(../../dlib/)

set(HDF5_USE_STATIC_LIBRARIES   ON)

include_directories(${HDF5_INCLUDE_DIRS})
add_library(ini ini.c INIReader.cpp)
add_library(Points Points.cpp)
add_library(utils utils.cpp)
add_library(tenXReader tenXReader.cpp)
#