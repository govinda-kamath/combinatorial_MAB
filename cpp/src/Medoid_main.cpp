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


int main(int argc, char *argv[]){
    std::string nameConfig = argv[1];

//    std::string nameConfig = "nominal.ini"; //For debugging mode in CLion

    INIReader reader(nameConfig);
    if (reader.ParseError() < 0) {
        std::cout << "Can't load "<< nameConfig << std::endl;
        return 1;
    }
    std::string directoryPath = reader.Get("path", "directory", "");
    std::string saveFilePath =reader.Get("path", "saveFilePath", "test.output");
    std::string fileSuffix = reader.Get("path", "suffix", "");
    unsigned numberOfInitialPulls = (unsigned) reader.GetInteger("UCB", "numberOfInitialPulls_medoid", 100);
    float delta = (float) reader.GetReal("UCB", "delta", 0.01);

    std::cout << "Running Medoid" << std::endl;

    // Load data
    std::vector<std::string>  pathsToImages;
    utils::getPathToFile(pathsToImages, directoryPath, fileSuffix);

    // Points
    std::vector<L1Point> pointsVec;
    utils::vectorsToPoints(pointsVec, pathsToImages);

    //Arms
    std::vector<ArmMedoid<L1Point> > armsVec;
    for (unsigned i(0); i < pointsVec.size(); i++) {
        ArmMedoid<L1Point> tmpArm(i, pointsVec[i], pointsVec);
        armsVec.push_back(tmpArm);
    }
    std::chrono::system_clock::time_point loopTimeStart = std::chrono::system_clock::now();
    std::cout<< "Starting UCB" <<std::endl;

    //Running UCB
    unsigned k = 10;
    UCB<ArmMedoid<L1Point> > UCB1(armsVec, delta, 3);
    UCB1.initialise(numberOfInitialPulls);
    UCB1.runUCB(20000*pointsVec.size());

    // True Means
    std::vector<float> topKArmsTrueMean(k);
    for (unsigned i = 0; i < k*5; i++) {
        topKArmsTrueMean[i] = UCB1.topKArms[i].trueMean();
    }
    std::vector<int> topKArmsArgSort(k);
    std::iota(topKArmsArgSort.begin(), topKArmsArgSort.end(), 0);
    auto comparator = [&topKArmsTrueMean](int a, int b){ return topKArmsTrueMean[a] < topKArmsTrueMean[b]; };
    std::sort(topKArmsArgSort.begin(), topKArmsArgSort.end(), comparator);


    //Print Result
    std::chrono::system_clock::time_point loopTimeEnd = std::chrono::system_clock::now();
    std::cout << "Average number of pulls " << UCB1.globalNumberOfPulls/UCB1.numberOfArms <<std::endl;
    std::cout << "Average time(ms) UCB "
    << std::chrono::duration_cast<std::chrono::milliseconds>(loopTimeEnd - loopTimeStart).count() << std::endl;
    std::cout<<"Rank\tId\tTr Mean\tEs Mean\tLCB\tUCB\tNoP"<<std::endl;
    for (unsigned i = 0; i < k; i++) {
        std::cout << topKArmsArgSort[i]+1
                << "\t" << UCB1.topKArms[i].id
                << "\t" << topKArmsTrueMean[i]
                << "\t" << UCB1.topKArms[i].estimateOfMean
                << "\t" << UCB1.topKArms[i].lowerConfidenceBound
                << "\t" << UCB1.topKArms[i].upperConfidenceBound
                << "\t" << UCB1.topKArms[i].numberOfPulls << std::endl;
    }

}
