cmake_minimum_required(VERSION 3.0.0)
project(combinatorial_MAB)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

option(knn "KNN" OFF)
option(medoid "Medoid" OFF)
option(kmeans "Kmeans" OFF)
option(med10x "Medoid10x" OFF)
option(knn10x "Knn10x" OFF)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(CMAKE_CXX_STANDARD 11)

#Bring the headers, such as Student.h into the project
include_directories(cpp/include)
include_directories(cpp/utilities)

#Can manually add the sources using the set command as follows:
#set(SOURCES src/mainapp.cpp src/Student.cpp)
add_subdirectory(dlib/dlib dlib_build)
add_subdirectory(nanoflann nanoflann )
add_subdirectory(cpp/lib lib)

find_package(HDF5 COMPONENTS C CXX HL REQUIRED)
#link_directories( ${HDF5_LIBRARY_DIRS} )
#include_directories( ${HDF5_INCLUDE_DIRS} )


file(GLOB INCLUDE "cpp/include/*.h")
file(GLOB LIBRARIES "cpp/lib/*.cpp")

set(Boost_USE_STATIC_LIBS   ON)

FIND_PACKAGE(Boost)
#FIND_PACKAGE ( Boost COMPONENTS container REQUIRED)
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )



#find_package(HDF5 REQUIRED COMPONENTS C CXX)
#add_executable(readdata cpp/deprecated/readdata.cpp cpp/include/tenXReader.h cpp/lib/tenXReader.cpp)
#target_link_libraries( readdata hdf5 hdf5_cpp)


#However, the file(GLOB...) allows for wildcard additions:
#file(GLOB SOURCES "cpp/src/*.cpp")
#add_executable(knn ${SOURCES} ${INCLUDE} ${LIBRARIES})

#if(${medoid})
#    add_executable(medoid cpp/src/Medoid_main.cpp
#    cpp/include/Points.h cpp/include/Arms.h
#    cpp/include/ini.h cpp/include/INIReader.h cpp/lib/ini.c
#    cpp/lib/Points.cpp
#    cpp/lib/INIReader.cpp cpp/include/utils.h cpp/lib/utils.cpp)
#    target_link_libraries(medoid dlib::dlib pthread)
#elseif(${knn})
#    add_executable(knn cpp/src/Knn_main.cpp
#        cpp/include/Points.h cpp/include/Arms.h
#        cpp/include/ini.h cpp/include/INIReader.h cpp/lib/ini.c
#        cpp/lib/Points.cpp
#        cpp/lib/INIReader.cpp cpp/include/utils.h cpp/lib/utils.cpp cpp/include/Knn.h
#            cpp/include/kmeans.h)
#    target_link_libraries(knn dlib::dlib pthread)
#elseif(${kmeans})
#    add_executable(kmeans cpp/src/Kmeans_main.cpp
#            cpp/include/Points.h cpp/include/Arms.h
#            cpp/include/ini.h cpp/include/INIReader.h cpp/lib/ini.c
#            cpp/lib/Points.cpp
#            cpp/lib/INIReader.cpp cpp/include/utils.h cpp/lib/utils.cpp
#            cpp/include/kmeans.h)
#    target_link_libraries(kmeans dlib::dlib pthread)
#elseif(${med10x})
#    add_executable(medoid10x cpp/src/Medoid_10x.cpp
#            cpp/include/Points.h cpp/include/Arms.h
#            cpp/include/ini.h cpp/include/INIReader.h cpp/lib/ini.c
#            cpp/lib/Points.cpp
#            cpp/lib/INIReader.cpp cpp/include/utils.h cpp/lib/utils.cpp
#            cpp/include/tenXReader.h cpp/lib/tenXReader.cpp)
#    target_link_libraries( medoid10x dlib::dlib ${HDF5_CXX_LIBRARIES})
#elseif(${knn10x})
#    add_executable(knn10x cpp/src/Knn_main_10x.cpp
#            cpp/include/Points.h cpp/include/Arms.h
#            cpp/include/ini.h cpp/include/INIReader.h cpp/lib/ini.c
#            cpp/lib/Points.cpp
#            cpp/lib/INIReader.cpp cpp/include/utils.h cpp/lib/utils.cpp
#            cpp/include/tenXReader.h cpp/lib/tenXReader.cpp)
#    target_link_libraries( knn10x dlib::dlib ${HDF5_CXX_LIBRARIES})

#else()
#    add_executable(knn cpp/src/Knn_Imagenet.cpp
#            cpp/include/Points.h cpp/include/Arms.h
#            cpp/utilities/ini.h cpp/utilities/INIReader.h cpp/utilities/ini.c
#            cpp/lib/Points.cpp
#            cpp/utilities/INIReader.cpp cpp/include/utils.h cpp/lib/utils.cpp cpp/include/Knn.h
#            cpp/include/kmeans.h )
#    target_link_libraries(knn dlib::dlib pthread ${Boost_LIBRARIES})
#
#    add_executable(kmeans cpp/src/Kmeans_Imagenet.cpp
#            cpp/include/Points.h cpp/include/Arms.h
#            cpp/utilities/ini.h cpp/utilities/INIReader.h cpp/utilities/ini.c
#            cpp/lib/Points.cpp
#            cpp/utilities/INIReader.cpp cpp/include/utils.h cpp/lib/utils.cpp
#            cpp/include/kmeans.h )
#    target_link_libraries(kmeans dlib::dlib pthread ${Boost_LIBRARIES})

    add_executable(heirarchical cpp/src/Heirarchical_Imagenet.cpp
            cpp/include/UCB_dynamic.h cpp/include/Points.h cpp/include/Arms.h
            cpp/utilities/ini.h cpp/utilities/INIReader.h cpp/utilities/ini.c
            cpp/lib/Points.cpp
            cpp/utilities/INIReader.cpp cpp/include/utils.h cpp/lib/utils.cpp
            )
    target_link_libraries( heirarchical dlib::dlib )

#    add_executable(medoid cpp/src/Medoid_Imagenet.cpp
#            cpp/include/Points.h cpp/include/Arms.h
#            cpp/utilities/ini.h cpp/utilities/INIReader.h cpp/utilities/ini.c
#            cpp/lib/Points.cpp
#            cpp/utilities/INIReader.cpp cpp/include/utils.h cpp/lib/utils.cpp)
#    target_link_libraries(medoid dlib::dlib pthread ${Boost_LIBRARIES})
#

#    add_executable(kmeans10x cpp/src/Kmeans_10x.cpp
#            cpp/include/Points.h cpp/include/Arms.h
#            cpp/utilities/ini.h cpp/utilities/INIReader.h cpp/utilities/ini.c
#            cpp/lib/Points.cpp
#            cpp/utilities/INIReader.cpp cpp/include/utils.h cpp/lib/utils.cpp
#            cpp/include/tenXReader.h cpp/lib/tenXReader.cpp)
#    target_link_libraries( kmeans10x dlib::dlib ${HDF5_CXX_LIBRARIES} ${Boost_LIBRARIES})
#
#
#    add_executable(knn10x cpp/src/Knn_10x.cpp
#            cpp/include/Points.h cpp/include/Arms.h
#            cpp/utilities/ini.h cpp/utilities/INIReader.h cpp/utilities/ini.c
#            cpp/lib/Points.cpp
#            cpp/utilities/INIReader.cpp cpp/include/utils.h cpp/lib/utils.cpp
#            cpp/include/tenXReader.h cpp/lib/tenXReader.cpp)
#    target_link_libraries( knn10x dlib::dlib ${HDF5_CXX_LIBRARIES} ${Boost_LIBRARIES})
#
#
##
    add_executable(heirarchical10x cpp/src/Heirarchical_10x.cpp
            cpp/include/UCB_dynamic.h cpp/include/Points.h cpp/include/Arms.h
            cpp/utilities/ini.h cpp/utilities/INIReader.h cpp/utilities/ini.c
            cpp/lib/Points.cpp cpp/include/tenXReader.h
            cpp/utilities/INIReader.cpp cpp/include/utils.h cpp/lib/utils.cpp
            cpp/lib/tenXReader.cpp
            )
    target_link_libraries( heirarchical10x dlib::dlib ${HDF5_CXX_LIBRARIES} ${Boost_LIBRARIES})

    add_executable(nndescent cpp/nndescent_wei/test_case.cpp
            cpp/nndescent_wei/nndes.h cpp/nndescent_wei/nndes-common.h
            cpp/include/UCB_dynamic.h cpp/include/Points.h cpp/include/Arms.h
            cpp/utilities/ini.h cpp/utilities/INIReader.h cpp/utilities/ini.c
            cpp/lib/Points.cpp cpp/include/tenXReader.h
            cpp/utilities/INIReader.cpp cpp/include/utils.h cpp/lib/utils.cpp
            cpp/lib/tenXReader.cpp
            )
    target_link_libraries(nndescent dlib::dlib pthread ${Boost_LIBRARIES} ${HDF5_CXX_LIBRARIES} nanoflann)

#    add_executable(medoid10x cpp/src/Medoid_10x.cpp
#            cpp/include/Points.h cpp/include/Arms.h
#            cpp/utilities/ini.h cpp/utilities/INIReader.h cpp/utilities/ini.c
#            cpp/lib/Points.cpp
#            cpp/utilities/INIReader.cpp cpp/include/utils.h cpp/lib/utils.cpp
#            cpp/include/tenXReader.h cpp/lib/tenXReader.cpp)
#    target_link_libraries( medoid10x dlib::dlib ${HDF5_CXX_LIBRARIES} ${Boost_LIBRARIES})

#    add_executable(knnGlove cpp/src/Knn_Glove.cpp
#            cpp/include/Points.h cpp/include/Arms.h
#            cpp/utilities/ini.h cpp/utilities/INIReader.h cpp/utilities/ini.c
#            cpp/lib/Points.cpp
#            cpp/utilities/INIReader.cpp cpp/include/utils.h cpp/lib/utils.cpp cpp/include/Knn.h
#            cpp/include/kmeans.h )
#    target_link_libraries(knnGlove dlib::dlib pthread ${Boost_LIBRARIES})

#endif()








