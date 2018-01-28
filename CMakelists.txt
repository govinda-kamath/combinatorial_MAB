cmake_minimum_required(VERSION 3.0.0)
project(combinatorial_MAB)

#Bring the headers, such as Student.h into the project
include_directories(cpp/include)

#Can manually add the sources using the set command as follows:
#set(SOURCES src/mainapp.cpp src/Student.cpp)

#However, the file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "cpp/src/*.cpp")

add_executable(test_code ${SOURCES})