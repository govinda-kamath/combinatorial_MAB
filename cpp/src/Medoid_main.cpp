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

    std::vector<std::string>  pathsToImages;
    utils::getPathToFile(pathsToImages, directoryPath, fileSuffix);

    std::vector<L1Point> pointsVec;
    std::vector<ArmMedoid<L1Point> > armsVec;

    // Points
    for  (unsigned long i(0); i < pathsToImages.size(); i++) {
        std::vector<float> tmpVec;
        utils::readImageAsVector(pathsToImages[i],tmpVec);
        L1Point tmpPoint(tmpVec);
        pointsVec.push_back(tmpPoint);
        if (i%10000 == 9999){
            std::cout << i+1 << " points read." << std::endl;
        }
    }
    std::cout<<"Data loaded" <<std::endl;
    //Arms
    for (unsigned i(0); i < pointsVec.size(); i++) {
        ArmMedoid<L1Point> tmpArm(i, pointsVec[i], pointsVec);
        armsVec.push_back(tmpArm);
    }




    std::chrono::system_clock::time_point loopTimeStart = std::chrono::system_clock::now();
    std::cout<< "Starting UCB" <<std::endl;

    //Running UCB
    unsigned k = 10;
    UCB<ArmMedoid<L1Point> > UCB1(armsVec, delta, 2);
    UCB1.initialise(numberOfInitialPulls);
    UCB1.runUCB(20000*pointsVec.size());

    //Result
    std::chrono::system_clock::time_point loopTimeEnd = std::chrono::system_clock::now();
    std::cout << "Average number of pulls " << UCB1.globalNumberOfPulls/UCB1.numberOfArms <<std::endl;
    std::cout << "Average time(ms) UCB "
    << std::chrono::duration_cast<std::chrono::milliseconds>(loopTimeEnd - loopTimeStart).count() << std::endl;
    std::cout<<"Id\tEsti\tLCB\tUCB\tNoP"<<std::endl;
    for (unsigned i = 0; i < k; i++) {
        std::cout << UCB1.topKArms[i].id
                << "\t" << UCB1.topKArms[i].estimateOfMean
                << "\t" << UCB1.topKArms[i].lowerConfidenceBound
                << "\t" << UCB1.topKArms[i].upperConfidenceBound
                << "\t" << UCB1.topKArms[i].numberOfPulls
                  << std::endl;
    }
    std::cout<<std::endl;

//    std::cout<<"The top estimates are";
//    for (unsigned i = 0; i < k; i++) {
//        std::cout <<  ;
//    }
//    std::cout<<std::endl;
//    std::cout<<"The top LCB are";
//    for (unsigned i = 0; i < k; i++) {
//        std::cout <<  " " << UCB1.topKArms[i].lowerConfidenceBound;
//    }
//    std::cout<<std::endl;


    // Naive
    std::vector<float> trueMean;
    for (unsigned i(0); i < pointsVec.size(); i++) {
        trueMean.push_back(armsVec[i].trueMean());
    }

    loopTimeEnd = std::chrono::system_clock::now();
    std::cout << "Average time(ms) Naive "
              << std::chrono::duration_cast<std::chrono::milliseconds>(loopTimeEnd - loopTimeStart).count() << std::endl;

    //Stats
    std::vector<int> trueMeanArgSort(pointsVec.size());
    std::iota(trueMeanArgSort.begin(), trueMeanArgSort.end(), 0);
    auto comparator = [&trueMean](int a, int b){ return trueMean[a] < trueMean[b]; };
    std::sort(trueMeanArgSort.begin(), trueMeanArgSort.end(), comparator);

    std::cout<< "Starting Naive method" <<std::endl;
    loopTimeStart = std::chrono::system_clock::now();

    std::cout<<"Top true Medoids are";
    for (unsigned i = 0; i < k; i++) {
        std::cout <<  " " << trueMeanArgSort[i];
    }
    std::cout<<std::endl;



    std::cout<<"The true means   are";
    for (unsigned i = 0; i < k; i++) {
        std::cout <<  " " << trueMean[trueMeanArgSort[i]];
    }
    std::cout<<std::endl;


    std::cout<<std::endl;
    return 0;
}
