cmake_minimum_required(VERSION 3.0)

include_directories(../../dlib/)

set(HDF5_USE_STATIC_LIBRARIES   ON)
add_library(ini ini.c INIReader.cpp)