//
// Created by Govinda Kamath on 2/8/18.
//

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <thread>
#include <numeric>
#include <mutex>
#include <map>
#include <chrono>
#include <dlib/image_io.h>
#include <dlib/image_transforms.h>
#include "tenXReader.h"
#include "Points.h"
#include "Arms.h"
#include "UCB.h"
#include "INIReader.h"
#include "utils.h"
#include "INIReader.h"

int main(int argc, char *argv[]){


    std::string nameConfig = argv[1];

    INIReader reader(nameConfig);
    if (reader.ParseError() < 0) {
        std::cout << "Can't load "<< nameConfig << std::endl;
        return 1;
    }

    // Reading parameters
    std::string saveFilePath =reader.Get("path", "saveFilePath", "test.output");
    unsigned numberOfInitialPulls = (unsigned) reader.GetInteger("UCB", "numberOfInitialPulls_medoid", 100);
    float delta = (float) reader.GetReal("UCB", "delta", 0.01);
    std::string fileName = reader.Get("path", "h5path", "test_dataset/1M_neurons_neuron20k.h5");
    unsigned k = 5;

    // Loading 10x data
    std::vector<unsigned long> shapeData(2);
    tenXReader::get10xMatrixSize(fileName, shapeData);
    std::vector<std::vector<float> > denseDataMatrix(shapeData[1], std::vector<float>(shapeData[0]));
    std::cout << "The dense matrix is " << denseDataMatrix.size() << " x " << denseDataMatrix[0].size() << std::endl;
    tenXReader::get10xNormalisedDenseMatrix(fileName, denseDataMatrix);

    //Arms
    std::vector<L1Point> pointsVec;
    utils::vectorsToPoints(pointsVec, denseDataMatrix);
    std::vector<ArmMedoid<L1Point> > armsVec;
    for (unsigned i(0); i < 100; i++) {
        ArmMedoid<L1Point> tmpArm(i, pointsVec[i], pointsVec);
        armsVec.push_back(tmpArm);
    }

    //UCB
    std::chrono::system_clock::time_point loopTimeStart = std::chrono::system_clock::now();
    std::cout<< "Starting UCB" <<std::endl;
    std::srand(std::time(nullptr));
    UCB<ArmMedoid<L1Point> > UCB1(armsVec, delta, k);

    std::cout << "best arm = " << UCB1.bestArm().id <<std::endl;
    // UCB: Initialization
    UCB1.initialise(numberOfInitialPulls);
    std::cout<< "Finished initialising" <<std::endl;

    for(int i(0); i < 100 ; i++){
        std::cout << i << " " << UCB1.armsContainer[i].lowerConfidenceBound
                  << " " << UCB1.armsContainer[i].estimateOfMean
                  << " " << UCB1.armsContainer[i].upperConfidenceBound
                  << " " << UCB1.armsContainer[i].numberOfPulls
                  << " " << UCB1.armsContainer[i].trueMean()
                  << std::endl;
    }
    loopTimeStart = std::chrono::system_clock::now();

    //
    std::cout << "Iterating" << std::endl;
    UCB1.runUCB(2000000*pointsVec.size());
    std::cout << "Iterating done" << std::endl;
    //Print Result
    std::chrono::system_clock::time_point loopTimeEnd = std::chrono::system_clock::now();
    std::cout << "Average number of pulls " << UCB1.globalNumberOfPulls/UCB1.numberOfArms <<std::endl;
    std::cout << "Average time(ms) UCB "
              << std::chrono::duration_cast<std::chrono::milliseconds>(loopTimeEnd - loopTimeStart).count() << std::endl;
    std::cout<<"Rank\tId\tTr Mean\tEs Mean\tLCB\tUCB\tNoP"<<std::endl;

    // True Means
    std::vector<float> topKArmsTrueMean(k*5);
    std::vector<int> topKArmsArgSort(k*5);
    for (unsigned i = 0; i < k*5; i++) {
        topKArmsTrueMean[i] = UCB1.topKArms[i].trueMean();
    }
    std::iota(topKArmsArgSort.begin(), topKArmsArgSort.end(), 0);
    auto comparator = [&topKArmsTrueMean](int a, int b){ return topKArmsTrueMean[a] < topKArmsTrueMean[b]; };
    std::sort(topKArmsArgSort.begin(), topKArmsArgSort.end(), comparator);


    for (unsigned i = 0; i < k*2; i++) {
        std::cout << std::setprecision (15) << topKArmsArgSort[i]+1
                  << "\t" << UCB1.topKArms[i].id
                  << "\t" << topKArmsTrueMean[i]
                  << "\t" << UCB1.topKArms[i].estimateOfMean
                  << "\t" << UCB1.topKArms[i].lowerConfidenceBound
                  << "\t" << UCB1.topKArms[i].upperConfidenceBound
                  << "\t" << UCB1.topKArms[i].numberOfPulls << std::endl;
    }

//    for(int i(7370); i < 7370+10 ; i++){
//        std::cout << i << " " << UCB1.armsContainer[i].lowerConfidenceBound
//                  << " " << UCB1.armsContainer[i].estimateOfMean
//                  << " " << UCB1.armsContainer[i].upperConfidenceBound
//                  << " " << UCB1.armsContainer[i].numberOfPulls
//                  << " " << UCB1.armsContainer[i].trueMean()
//                  << std::endl;
//    }
}
