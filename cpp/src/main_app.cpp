#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <thread>
#include <map>
#include <dlib/image_io.h>
#include <dlib/image_transforms.h>
#include "Points.h"
#include "Arms.h"
#include "UCB.h"
#include "INIReader.h"
#include "utils.h"
#define DEBUG


void singleRun(std::vector<SquaredEuclideanPoint> &pointsVec, unsigned threadNumber, unsigned long mainPointIndexStart,
               unsigned long mainPointIndexEnd, int numberOfInitialPulls, float delta){

    std::cout << "Running thread "<< threadNumber << std::endl;
    for (unsigned long index = mainPointIndexStart; index<mainPointIndexEnd; index++){
        std::chrono::system_clock::time_point loopTimeStart = std::chrono::system_clock::now();
//        std::this_thread::sleep_for(std::chrono::seconds(1*(threadNumber+1)));


        std::vector<ArmKNN<SquaredEuclideanPoint> > armsVec;
        for (unsigned i(0); i < pointsVec.size(); i++) {
            if (i == index)
                continue;
            ArmKNN<SquaredEuclideanPoint> tmpArm(i, pointsVec[i], pointsVec[index]);
            armsVec.push_back(tmpArm);
        }
        std::chrono::system_clock::time_point loopTimeMid = std::chrono::system_clock::now();

        UCB<ArmKNN<SquaredEuclideanPoint> > UCB1(armsVec, delta);
        UCB1.initialise(numberOfInitialPulls);
        UCB1.runUCB(1000000);


        std::chrono::system_clock::time_point loopTimeEnd = std::chrono::system_clock::now();
        if (index%5==0) {
            std::cout << "Ran thread " << threadNumber << ".Index " << index << " Total Time taken = " <<
                      std::chrono::duration_cast<std::chrono::milliseconds>(loopTimeEnd - loopTimeStart).count()
                      << std::endl;
            std::cout << "Ran thread " << threadNumber << ".Index " << index << " Mid Time taken = " <<
                      std::chrono::duration_cast<std::chrono::milliseconds>(loopTimeMid - loopTimeStart).count()
                      << std::endl;
        }
    }
}



int main(int argc, char *argv[]){
    std::cout << "We have entered " << argc
              << " arguments." << std::endl;

    std::vector<SquaredEuclideanPoint> pointsVec;

    std::string nameConfig;
    if (argc==2){
        std::cout<<argc<<std::endl;
        nameConfig = argv[1];
    }
    else{
        nameConfig = "nominal.ini";
    }
    INIReader reader(nameConfig);

    if (reader.ParseError() < 0) {
        std::cout << "Can't load "<< nameConfig << std::endl;
        return 1;
    }

    std::string directoryPath = reader.Get("path", "directory", "");
    std::string fileSuffix = reader.Get("path", "suffix", "");

    int numberOfInitialPulls = (int) reader.GetInteger("UCB", "numberOfInitialPulls", 100);
    float delta = (float) reader.GetReal("UCB", "delta", 0.1);
    int numCores = (int) reader.GetReal("UCB", "numCores", 1);
    unsigned long  maxNumberOfPoints = (unsigned long) reader.GetReal("UCB", "noPoints", 1);

    std::cout << "Running K-nn for " << maxNumberOfPoints << " points using ";
    std::cout << "Number of cores = " << numCores<<std::endl;
    std::cout << numberOfInitialPulls << std::endl;
    std::cout << delta << std::endl;
    std::cout << directoryPath << std::endl;
    std::cout << fileSuffix << std::endl;
    std::vector<float> tmpVec;

    unsigned long fileNumber(0);
    std::vector<std::string>  pathsToImages;
    clock_t timeRead = clock();
    utils::getPathToFile(pathsToImages, directoryPath, fileSuffix);

    unsigned long pointIndex(0);

    for  (unsigned long i(0); i < pathsToImages.size(); i++) {
        float tmpValue;

        std::vector<float> tmpVec;
        utils::readImageAsVector(pathsToImages[i],tmpVec);

        SquaredEuclideanPoint tmpPoint(tmpVec);
        pointsVec.push_back(tmpPoint);
        pointIndex++;

        if (pointIndex%2000 == 999){
            std::cout << pointIndex+1 << " points read." << std::endl;
        }
    }

    std::cout << "Reading time (ms)" << 1000 * (clock() - timeRead) / CLOCKS_PER_SEC << std::endl;

    //Parallelize
    std::chrono::system_clock::time_point loopTimeStart = std::chrono::system_clock::now();
    std::vector<std::thread> initThreads(numCores);
    unsigned long chunkSize = (maxNumberOfPoints/numCores);

//    #pragma omp parallel for
    for(unsigned t = 0; t < numCores; t++) {

        unsigned long mainPointIndexStart = t * chunkSize;
        unsigned long mainPointIndexEnd = (t + 1) * chunkSize;
        initThreads[t] = std::thread(singleRun, std::ref(pointsVec), t, mainPointIndexStart, mainPointIndexEnd, numberOfInitialPulls, delta);
//        singleRun(pointsVec, mainPointIndexStart, mainPointIndexEnd, numberOfInitialPulls, delta);


    }
    for(unsigned t = 0; t < numCores; t++) {
        initThreads[t].join();
    }
    std::chrono::system_clock::time_point loopTimeEnd = std::chrono::system_clock::now();

    std::cout << "Average time (ms)"
              << std::chrono::duration_cast<std::chrono::milliseconds>(loopTimeEnd - loopTimeStart).count()/maxNumberOfPoints << std::endl;

    return 0;
}