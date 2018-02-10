//
// Created by Govinda Kamath on 2/8/18.
//

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <thread>
#include <numeric>
#include <mutex>
#include <map>
#include<chrono>
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

    std::string saveFilePath =reader.Get("path", "saveFilePath", "test.output");
    unsigned numberOfInitialPulls = (unsigned) reader.GetInteger("UCB", "numberOfInitialPulls_medoid", 100);
    float delta = (float) reader.GetReal("UCB", "delta", 0.01);

    std::string fileName = reader.Get("path", "h5path", "test_dataset/1M_neurons_neuron20k.h5");

    std::vector<int> shapeData(2);
    tenXReader::get10xMatrixSize(fileName, shapeData);

    std::chrono::system_clock::time_point dataLoadTimeStart = std::chrono::system_clock::now();
    std::vector<std::vector<float> > denseDataMatrix(shapeData[1], std::vector<float>(shapeData[0]));
    std::cout << denseDataMatrix.size() << std::endl;
    std::cout << denseDataMatrix[0].size() << std::endl;
    std::cout << denseDataMatrix[0][100] << std::endl;

    tenXReader::get10xNormalisedDenseMatrix(fileName, denseDataMatrix);
    std::chrono::system_clock::time_point dataLoadTimeEnd = std::chrono::system_clock::now();
    std::cout << "Time to read dense vector (ms) "
              << std::chrono::duration_cast<std::chrono::milliseconds>(dataLoadTimeEnd - dataLoadTimeStart).count()
              << std::endl;

    for (int i = 0; i < 10; ++i) {
        int numberNonzeros(0);
        for (int j = 0; j < denseDataMatrix[0].size(); ++j){
            if (denseDataMatrix[i][j] != 0){
                numberNonzeros++;
            }
        }
        std::cout << i << "  " << numberNonzeros << std::endl;
    }

    dataLoadTimeStart = std::chrono::system_clock::now();
//    std::cout << "row 0 size " << denseDataMatrix[0].size() << std::endl;
//    std::cout << "row 120 size " << denseDataMatrix[120].size() << std::endl;
//    std::cout << "row number  " << denseDataMatrix.size() << std::endl;
    std::vector<L1Point> pointsVec;
    utils::vectorsToPoints(pointsVec, denseDataMatrix);

//    std::cout << "arm 0 size " << pointsVec[0].vecSize << std::endl;
//    std::cout << "row 120 size " << pointsVec[120].vecSize << std::endl;
//    std::cout << "row number  " << pointsVec.size() << std::endl;

    dataLoadTimeEnd = std::chrono::system_clock::now();
    std::cout << "Time to convert to points (ms) "
              << std::chrono::duration_cast<std::chrono::milliseconds>(dataLoadTimeEnd - dataLoadTimeStart).count()
              << std::endl;

    //Arms
    std::vector<ArmMedoid<L1Point> > armsVec;
    for (unsigned i(0); i < pointsVec.size(); i++) {
        ArmMedoid<L1Point> tmpArm(i, pointsVec[i], pointsVec);
//        std::cout << "i = " << i << " size " << pointsVec[i].point.size() << std::endl ;
        armsVec.push_back(tmpArm);
    }
    std::chrono::system_clock::time_point loopTimeStart = std::chrono::system_clock::now();
    std::cout<< "Starting UCB" <<std::endl;

    //Running UCB
    unsigned k = 5;
    UCB<ArmMedoid<L1Point> > UCB1(armsVec, delta, k);
    std::cout << numberOfInitialPulls << "pulls initially" << std::endl;
    UCB1.initialise(numberOfInitialPulls);
    std::cout<< "Finished initialising" <<std::endl;

    for(int i(7370); i < 7370+10 ; i++){
        std::cout << i << " " << UCB1.armsContainer[i].lowerConfidenceBound
                  << " " << UCB1.armsContainer[i].estimateOfMean
                  << " " << UCB1.armsContainer[i].upperConfidenceBound
                  << " " << UCB1.armsContainer[i].numberOfPulls
                  << " " << UCB1.armsContainer[i].trueMean()
                  << std::endl;
    }
    loopTimeStart = std::chrono::system_clock::now();
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

    for(int i(7370); i < 7370+10 ; i++){
        std::cout << i << " " << UCB1.armsContainer[i].lowerConfidenceBound
                  << " " << UCB1.armsContainer[i].estimateOfMean
                  << " " << UCB1.armsContainer[i].upperConfidenceBound
                  << " " << UCB1.armsContainer[i].numberOfPulls
                  << " " << UCB1.armsContainer[i].trueMean()
                  << std::endl;
    }
}
