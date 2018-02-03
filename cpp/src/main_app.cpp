#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <thread>
#include <numeric>
#include <map>
#include <dlib/image_io.h>
#include <dlib/image_transforms.h>
#include "Points.h"
#include "Arms.h"
#include "UCB.h"
#include "INIReader.h"
#include "utils.h"
#define DEBUG



void singleRun(std::vector<SquaredEuclideanPoint> &pointsVec, unsigned  k, unsigned long mainPointIndexStart,
               unsigned long mainPointIndexEnd, unsigned numberOfInitialPulls, float delta, std::string saveFilepath,
               std::vector<std::string>  pathsToImages){


    std::ofstream saveFile;
    saveFile.open (saveFilepath, std::ofstream::out | std::ofstream::app);
    std::vector<std::vector<ArmKNN<SquaredEuclideanPoint>> > allAnswers;
    std::vector<float> avgNumberOfPulls;

    for (unsigned long index = mainPointIndexStart; index<mainPointIndexEnd; index++){
        std::vector<ArmKNN<SquaredEuclideanPoint> > armsVec;
        for (unsigned i(0); i < pointsVec.size(); i++) {
            if (i == index)
                continue;
            ArmKNN<SquaredEuclideanPoint> tmpArm(i, pointsVec[i], pointsVec[index]);
            armsVec.push_back(tmpArm);
        }

        UCB<ArmKNN<SquaredEuclideanPoint> > UCB1(armsVec, delta, k);
        UCB1.initialise(numberOfInitialPulls);
        UCB1.runUCB(200*pointsVec.size());
        allAnswers.push_back(UCB1.topKArms);
        if (index%100==0){
            std::cout << "Thread " << mainPointIndexStart << ". Index " << index<< " "
                      << pathsToImages[UCB1.topKArms[0].id]<<std::endl;
        }
        avgNumberOfPulls.push_back(UCB1.globalNumberOfPulls/UCB1.numberOfArms);
    }
    std::cout<< "Saving the thread starting with " << mainPointIndexStart << " in file " << saveFilepath << std::endl;

    for (unsigned long index = mainPointIndexStart; index<mainPointIndexEnd ; index++) {

        std::vector<ArmKNN<SquaredEuclideanPoint>> topKArms = allAnswers[index - mainPointIndexStart];
        std::vector<float> topKArmsTrueMean(k);

        saveFile << index << "R ";
        for (unsigned i = 0; i < k*5; i++) {
            saveFile << topKArms[i].id << " ";
            topKArmsTrueMean[i] = topKArms[i].trueMean();
        }
        saveFile << std::endl;

        std::vector<int> topKArmsArgSort(k*5);
        std::iota(topKArmsArgSort.begin(), topKArmsArgSort.end(), 0);
        auto comparator = [&topKArmsTrueMean](int a, int b){ return topKArmsTrueMean[a] < topKArmsTrueMean[b]; };
        std::sort(topKArmsArgSort.begin(), topKArmsArgSort.end(), comparator);

        saveFile << index << "A";
        for (unsigned i = 0; i < k*5; i++) {
            saveFile <<  " " << topKArmsArgSort[i];
        }
        saveFile << "Av:" << avgNumberOfPulls[index - mainPointIndexStart] << "\n";
        saveFile << std::endl;

    }

}



int main(int argc, char *argv[]){
    std::string nameConfig = argv[1];
    long startIndex(atol(argv[2])); // Start index
    long endIndex(atol(argv[3])); // End index


//     For debugging mode in CLion
//    std::string nameConfig = "nominal.ini";
//    long startIndex(0); // Start index
//    long endIndex(100); // End index


    INIReader reader(nameConfig);
    if (reader.ParseError() < 0) {
        std::cout << "Can't load "<< nameConfig << std::endl;
        return 1;
    }
    std::string directoryPath = reader.Get("path", "directory", "");
    std::string saveFilePath =reader.Get("path", "saveFilePath", "test.output");
    std::string fileSuffix = reader.Get("path", "suffix", "");
    unsigned numberOfInitialPulls = (unsigned) reader.GetInteger("UCB", "numberOfInitialPulls", 100);
    unsigned k = (unsigned) reader.GetInteger("UCB", "k", 5);
    float delta = (float) reader.GetReal("UCB", "delta", 0.1);
//    int numCores = (int) reader.GetReal("UCB", "numCores", 1);

    std::cout << "Running "<<k<<"nn for " << endIndex - startIndex << std::endl;
    std::cout << numberOfInitialPulls << std::endl;
    std::cout << delta << std::endl;
    std::cout << directoryPath << std::endl;
    std::vector<float> tmpVec;

    std::vector<std::string>  pathsToImages;
    clock_t timeRead = clock();
    utils::getPathToFile(pathsToImages, directoryPath, fileSuffix);

    unsigned long pointIndex(0);
    std::vector<SquaredEuclideanPoint> pointsVec;
    for  (unsigned long i(0); i < pathsToImages.size(); i++) {
        std::vector<float> tmpVec;
        utils::readImageAsVector(pathsToImages[i],tmpVec);
        SquaredEuclideanPoint tmpPoint(tmpVec);
        pointsVec.push_back(tmpPoint);
        pointIndex++;
        if (pointIndex%10000 == 9999){
            std::cout << pointIndex+1 << " points read." << std::endl;
        }
    }

    std::cout << "Reading time (ms)" << 1000 * (clock() - timeRead) / CLOCKS_PER_SEC << std::endl;

    std::chrono::system_clock::time_point loopTimeStart = std::chrono::system_clock::now();
    singleRun(pointsVec, k ,startIndex, endIndex, numberOfInitialPulls, delta, saveFilePath, pathsToImages);
//    //Parallelize
//    std::vector<std::thread> initThreads(numCores);
//    unsigned long chunkSize = (maxNumberOfPoints/numCores);
//
////    #pragma omp parallel for
//    for(unsigned t = 0; t < numCores; t++) {
//
//        unsigned long mainPointIndexStart = t * chunkSize;
//        unsigned long mainPointIndexEnd = (t + 1) * chunkSize;
//        initThreads[t] = std::thread(singleRun, std::ref(pointsVec), t, mainPointIndexStart, mainPointIndexEnd, numberOfInitialPulls, delta);
////
//
//
//    }
//    for(unsigned t = 0; t < numCores; t++) {
//        initThreads[t].join();
//    }
    std::chrono::system_clock::time_point loopTimeEnd = std::chrono::system_clock::now();
    std::cout << "Average time (ms) "
              << std::chrono::duration_cast<std::chrono::milliseconds>(loopTimeEnd - loopTimeStart).count()/
                      (endIndex-startIndex) << std::endl;

    return 0;
}
