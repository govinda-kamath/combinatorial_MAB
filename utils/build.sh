
#!/bin/bash
pwd=$PWD

mkdir -p build
cd $pwd/build
cmake .. -DCMAKE_C_COMPILER=gcc-4.8 -DCMAKE_CXX_COMPILER=g++-4.8
make 
#make install

exit $?
