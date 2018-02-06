cmake_minimum_required(VERSION 3.0.0)
project(combinatorial_MAB)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

option(knn "KNN" ON)
option(medoid "Medoid" ON)
option(kmeans "Medoid" ON)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

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
#file(GLOB SOURCES "cpp/src/*.cpp")
#add_executable(knn ${SOURCES} ${INCLUDE} ${LIBRARIES})
if(${knn})
    add_executable(knn cpp/src/Knn_main.cpp
        cpp/include/Points.h cpp/include/Arms.h
        cpp/include/ini.h cpp/include/INIReader.h cpp/lib/ini.c
        cpp/lib/Points.cpp
        cpp/lib/INIReader.cpp cpp/include/utils.h cpp/lib/utils.cpp cpp/include/Knn.h
            cpp/include/kmeans.h)
    target_link_libraries(knn dlib::dlib pthread)
endif()

if(${medoid})
    add_executable(medoid cpp/src/Medoid_main.cpp
        cpp/include/Points.h cpp/include/Arms.h
        cpp/include/ini.h cpp/include/INIReader.h cpp/lib/ini.c
        cpp/lib/Points.cpp
        cpp/lib/INIReader.cpp cpp/include/utils.h cpp/lib/utils.cpp)
    target_link_libraries(medoid dlib::dlib pthread)
endif()

if(${kmeans})
    add_executable(kmeans cpp/src/Kmeans_main.cpp
            cpp/include/Points.h cpp/include/Arms.h
            cpp/include/ini.h cpp/include/INIReader.h cpp/lib/ini.c
            cpp/lib/Points.cpp
            cpp/lib/INIReader.cpp cpp/include/utils.h cpp/lib/utils.cpp
            cpp/include/kmeans.h)
    target_link_libraries(kmeans dlib::dlib pthread)
endif()