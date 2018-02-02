cmake_minimum_required(VERSION 3.0)

include_directories(../../dlib/)

add_library(ini ini.c INIReader.cpp)
add_library(Points Points.cpp)
add_library(utils utils.cpp)
#add_library(Arms Arms.cpp)
