
#!/bin/bash
pwd=$PWD

mkdir -p build
cd $pwd/build

if [ "$1" == "knn" ]
then
    cmake .. -DCMAKE_C_COMPILER=gcc-4.9 -DCMAKE_CXX_COMPILER=g++-4.9 -DMedoid=OFF
elif [ "$1" == "medoid" ]
then
    cmake .. -DCMAKE_C_COMPILER=gcc-4.9 -DCMAKE_CXX_COMPILER=g++-4.9 -DKNN=OFF
else
    cmake .. -DCMAKE_C_COMPILER=gcc-4.9 -DCMAKE_CXX_COMPILER=g++-4.9
fi
make 
#make install

exit $?
