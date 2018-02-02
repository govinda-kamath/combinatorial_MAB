cmake_minimum_required(VERSION 3.0.0)
project(combinatorial_MAB)
set(CMAKE_CXX_STANDARD 11)

#Bring the headers, such as Student.h into the project
include_directories(cpp/include)

#Can manually add the sources using the set command as follows:
#set(SOURCES src/mainapp.cpp src/Student.cpp)
add_subdirectory(dlib/dlib dlib_build)
add_subdirectory(cpp/lib lib)

file(GLOB INCLUDE "cpp/include/*.h")
file(GLOB LIBRARIES "cpp/lib/*.cpp")

#However, the file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "cpp/src/*.cpp")

#add_executable(test_code ${SOURCES} ${INCLUDE} ${LIBRARIES})
add_executable(test_code ${SOURCES} cpp/include/Points.h cpp/include/Arms.h
        cpp/include/ini.h cpp/include/INIReader.h cpp/lib/ini.c
        cpp/lib/Points.cpp
        cpp/lib/INIReader.cpp cpp/include/utils.h cpp/lib/utils.cpp)
#target_link_libraries(test_code cpp/include/Test)
target_link_libraries(test_code dlib::dlib pthread)
