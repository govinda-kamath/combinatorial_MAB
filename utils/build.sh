#!/bin/bash
pwd=$PWD

mkdir -p build
cd $pwd/build

if [ "$1" == "knn" ]
then
    cmake .. -DCMAKE_C_COMPILER=gcc-4.9 -DCMAKE_CXX_COMPILER=g++-4.9 -DKNN=ON
elif [ "$1" == "medoid" ]
then
    cmake .. -DCMAKE_C_COMPILER=gcc-4.9 -DCMAKE_CXX_COMPILER=g++-4.9 -DMedoid=ON
elif [ "$1" == "kmeans" ]
then
    cmake .. -DCMAKE_C_COMPILER=gcc-4.9 -DCMAKE_CXX_COMPILER=g++-4.9 -DKmeans=ON
elif [ "$1" == "medoid10x" ]
then
    cmake .. -DCMAKE_C_COMPILER=gcc-4.9 -DCMAKE_CXX_COMPILER=g++-4.9 -DMedoid10x=ON
elif [ "$1" == "knn10x" ]
then
    cmake .. -DCMAKE_C_COMPILER=gcc-4.9 -DCMAKE_CXX_COMPILER=g++-4.9 -DKnn10x=ON
else
    cmake .. -DCMAKE_C_COMPILER=gcc-4.9 -DCMAKE_CXX_COMPILER=g++-4.9
fi
make 
#make install

exit $?
