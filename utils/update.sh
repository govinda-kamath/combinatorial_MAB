#!/bin/bash

if [ "$2" == "scripts" ];
then rsync -rizP --delete --exclude '.*' --exclude '*.pyc' --exclude 'figures' python/ $1@shannon.stanford.edu:/data/MAB/work/code/python
fi

if [ "$2" == "utils" ];
then rsync -rizP --delete --exclude '.*' --exclude '*.pyc' --exclude 'figures' utils/ $1@shannon.stanford.edu:/data/MAB/work/code/utils
fi

if [ "$2" == "push" ];
then rsync -rizP --delete --exclude '.*' --exclude 'build' cpp/ $1@shannon.stanford.edu:/data/MAB/work/code/cpp/
fi

if [ "$2" == "pull" ];
then rsync -rizP --delete --exclude '.*' --exclude 'build' $1@shannon.stanford.edu:/data/MAB/work/code/cpp/ cpp/
fi

if [ "$2" == "makefile" ];
then rsync -rizP --delete --exclude '.*' --exclude 'build' CMakeLists.txt $1@shannon.stanford.edu:/data/MAB/work/code/CMakeLists.txt 
fi

if [ "$2" == "update" ];
then ssh -t $1@shannon.stanford.edu "export TEMP=/home/$1/tmp && cd /data/MAB/work/code/ && ./utils/build.sh"
fi

if [ "$2" == "all" ];
then rsync -rizP --delete --exclude '.*' --exclude 'data' --exclude '*.pyc' --exclude 'figures' --exclude 'build' --exclude 'cmake-build-debug' . $1@shannon.stanford.edu:/data/MAB/work/code/
fi