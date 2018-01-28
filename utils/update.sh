#!/bin/bash

if [ "$2" == "scripts" ];
then rsync -rizP --delete --exclude '.*' --exclude '*.pyc' --exclude 'figures' python/ $1@shannon.stanford.edu:/data/MAB/work/code/python
fi

if [ "$2" == "utils" ];
then rsync -rizP --delete --exclude '.*' --exclude '*.pyc' --exclude 'figures' utils/ $1@shannon.stanford.edu:/data/MAB/work/code/python/utils
fi

if [ "$2" == "push" ];
then rsync -rizP --delete --exclude '.*' --exclude 'build' cpp/ $1@shannon.stanford.edu:/home/$1/AwesomeAssembler/cpp/
fi

if [ "$2" == "pull" ];
then rsync -rizP --delete --exclude '.*' --exclude 'build' $1@shannon.stanford.edu:/data/MAB/work/code/cpp/ cpp/
fi

if [ "$2" == "all" ];
then rsync -rizP --delete --exclude '.*' --exclude 'data' --exclude '*.pyc' --exclude 'figures' --exclude 'build' . $1@shannon.stanford.edu:/data/MAB/work/code/
