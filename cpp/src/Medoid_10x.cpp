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
    std::srand(std::time(nullptr));

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
    unsigned k = 3;

    // Loading 10x data shape
    std::vector<unsigned> shapeData =  tenXReader::get10xMatrixSize(fileName);

    // Loading Dense matrix
//    std::vector<std::vector<float> > denseDataMatrix(shapeData[1], std::vector<float>(shapeData[0]));
//    std::cout << "The dense matrix is " << denseDataMatrix.size() << " x " << denseDataMatrix[0].size() << std::endl;
//    tenXReader::get10xNormalisedDenseMatrix(fileName, denseDataMatrix);

//     Loading Sparse matrix
    std::cout << "Reading normalised data sparsely " << std::endl;
    std::vector<std::unordered_map<unsigned long, float> > sparseNormalisedDataMatrix(shapeData[1] );
    tenXReader::get10xNormalisedSparseMatrix(fileName, sparseNormalisedDataMatrix);


    //Arms
//
    std::vector<SparseL1Point> pointsVec;
    utils::unorderedMapToPoints(pointsVec, sparseNormalisedDataMatrix, shapeData[0]);
    std::vector<ArmMedoid<SparseL1Point> > armsVec;

//    std::vector<L1Point> pointsVec;
//    utils::vectorsToPoints(pointsVec, denseDataMatrix);
//    std::vector<ArmMedoid<L1Point> > armsVec;

    unsigned numberOfPoints = pointsVec.size();
    for (unsigned i(0); i < numberOfPoints; i++) {
        ArmMedoid<SparseL1Point> tmpArm(i, pointsVec[i], pointsVec);
        armsVec.push_back(tmpArm);
    }

    //UCB
    std::cout<< "Starting UCB" <<std::endl;
    UCB<ArmMedoid<SparseL1Point> > UCB1(armsVec, delta, k);
    std::cout << "no of arms" << UCB1.armsContainer.size() << std::endl;

    std::cout << "arm point length" << UCB1.armsContainer[0].point->vecSize << std::endl;

//    std::cout << "best arm = " << UCB1.bestArm().id <<std::endl;
    // UCB: Initialization
    std::cout << "UCB: Initializing with " << numberOfInitialPulls << " points" << std::endl;
    std::chrono::system_clock::time_point loopTimeStart = std::chrono::system_clock::now();
    UCB1.initialise(numberOfInitialPulls);
    std::chrono::system_clock::time_point loopTimeEnd = std::chrono::system_clock::now();
    std::cout<< "UCB: Finished initialising" <<std::endl;
    std::cout << "Time taken (s) "
              << std::chrono::duration_cast<std::chrono::seconds>(loopTimeEnd - loopTimeStart).count() << std::endl;

//    std::vector<float> v(numberOfPoints);
//    float localvaravg(0);
//    for(int i(0); i < numberOfPoints ; i++){
//        v[i] = UCB1.armsContainer[i].trueMean();
//        localvaravg += UCB1.armsContainer[i].estimateOfSecondMoment -std::pow(UCB1.armsContainer[i].estimateOfMean,2);
//        std::cout << i << " " << UCB1.armsContainer[i].lowerConfidenceBound
//                  << " " << UCB1.armsContainer[i].estimateOfMean
//                  << " " << UCB1.armsContainer[i].upperConfidenceBound
//                  << " " << UCB1.armsContainer[i].numberOfPulls
//                  << " " << UCB1.armsContainer[i].trueMean()
//                  << std::endl;
//    }
//    double sum = std::accumulate(std::begin(v), std::end(v), 0.0);
//    double m =  sum / v.size();
//
//    double accum = 0.0;
//    std::for_each (std::begin(v), std::end(v), [&](const double d) {
//        accum += (d - m) * (d - m);
//    });
//
//    double var = accum / (v.size()-1);
//
//    std::cout<< "Global var " << UCB1.globalSigma*UCB1.globalSigma  << "\n"
//             << " Local var avg " << localvaravg/numberOfPoints << "\n"
//            << " Variance of mu "  << var
//            << " Mean of mu "  << m
//            << std::endl;

//    loopTimeStart = std::chrono::system_clock::now();

    //
    std::cout << "UCB: Iterating" << std::endl;
    loopTimeStart = std::chrono::system_clock::now();
    UCB1.runUCB(100000*pointsVec.size());
    std::cout << "UCB: Iterating done" << std::endl;
    loopTimeEnd = std::chrono::system_clock::now();
    std::cout << "Time taken (s) "
              << std::chrono::duration_cast<std::chrono::seconds>(loopTimeEnd - loopTimeStart).count() << std::endl;

    //Print Result
    std::cout << "Total number of pulls " << UCB1.globalNumberOfPulls <<std::endl;
    std::cout << "Number of points " << UCB1.numberOfArms << std::endl;
    std::cout << "Dimension of each point" << shapeData[0] << std::endl;
    std::cout << "Average " << UCB1.globalNumberOfPulls/(UCB1.numberOfArms*shapeData[0]*0.07) <<std::endl;//0.07=sparsity
    std::cout << "Global Sigma = " << UCB1.globalSigma <<std::endl;

    // True Means
    std::vector<float> topKArmsTrueMean(std::min(k*2, numberOfPoints));
    std::vector<int> topKArmsArgSort(std::min(k*2, numberOfPoints));
    for (unsigned i = 0; i < std::min(k*2, numberOfPoints); i++) {
        topKArmsTrueMean[i] = UCB1.topKArms[i].trueMean();
    }
    std::iota(topKArmsArgSort.begin(), topKArmsArgSort.end(), 0);
    auto comparator = [&topKArmsTrueMean](int a, int b){ return topKArmsTrueMean[a] < topKArmsTrueMean[b]; };
    std::sort(topKArmsArgSort.begin(), topKArmsArgSort.end(), comparator);

    std::cout<<"Rank\tId\tTr Mean\t\t\t\t\tEs Mean\t\t\t\tLCB\t\t\t\t\t\tUCB\t\t\t\t\t\tNoP"<<std::endl;
    for (unsigned i = 0; i < std::min(k*2, numberOfPoints); i++) {
        std::cout << std::setprecision (15) << topKArmsArgSort[i]+1
                  << "\t\t" << UCB1.topKArms[i].id
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
