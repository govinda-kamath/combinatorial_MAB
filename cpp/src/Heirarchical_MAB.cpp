//
// Created by Govinda Kamath on 2/28/18.
//

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <thread>
#include <numeric>
#include <map>
#include "Points.h"
#include "Arms.h"
#include "../utilities/INIReader.h"
#include "Knn.h"


int main() {

//    std::string nameConfig = argv[1];
//    long startIndex(atol(argv[2])); // Start index
//    long endIndex(atol(argv[3])); // End index

//     For debugging mode in CLion
    std::string nameConfig = "/Users/govinda/Code/combinatorial_MAB/nominal1.ini";
    long startIndex(0); // Start index
    long endIndex(10); // End index

    // Parameters
    INIReader reader(nameConfig);
    if (reader.ParseError() < 0) {
        std::cout << "Can't load " << nameConfig << std::endl;
        return 1;
    }
    std::string directoryPath = reader.Get("path", "directory", "");
    std::string saveFilePath = reader.Get("path", "saveFilePath", "test.output");
    std::string fileSuffix = reader.Get("path", "suffix", "");
    unsigned numberOfInitialPulls = (unsigned) reader.GetInteger("UCB", "numberOfInitialPulls_knn", 100);
    float delta = (float) reader.GetReal("UCB", "delta", 0.1);

    std::cout << "Running Hierar" << "-nn for " << endIndex - startIndex << " points" << std::endl;
    std::cout << numberOfInitialPulls << std::endl;

    std::vector<std::string> pathsToImages;
    std::vector<SquaredEuclideanPoint> pointsVec;
    std::vector<std::shared_ptr<SquaredEuclideanPoint>> sharedPtrPointsVec;
    std::vector<GroupPoint<SquaredEuclideanPoint> > groupPoints;
    std::vector<ArmHeirarchical<SquaredEuclideanPoint> > armsVec;
    // Stores map from group id to arm id. Used to removes "irrelevant" arms
    std::unordered_map<std::pair<unsigned long, unsigned long>, unsigned long, utils::pair_hash> groupIDtoArmID;


    utils::getPathToFile(pathsToImages, directoryPath, fileSuffix);
    utils::vectorsToPoints(pointsVec, pathsToImages);
    groupIDtoArmID.reserve(armsVec.size() * armsVec.size());

    // Step 1: Initialize all the leaves
    unsigned long armID = 0;
    unsigned long n = pointsVec.size();
    unsigned long d = pointsVec[0].getVecSize();
    unsigned long maxPointId;

    for (unsigned long i(0); i < n; i++) {
        sharedPtrPointsVec.push_back(std::make_shared<SquaredEuclideanPoint>(pointsVec[i]));
        std::vector<std::shared_ptr<SquaredEuclideanPoint> > groupPointTmp1;
        groupPointTmp1.push_back(sharedPtrPointsVec[i]);
        GroupPoint<SquaredEuclideanPoint> gpTmp1(groupPointTmp1, d, i);
        groupPoints.push_back(gpTmp1);
    }

    maxPointId = n-1;

    for (unsigned long i(0); i < n; i++) {
        for (unsigned long j(0); j < i; j++) {
            ArmHeirarchical<SquaredEuclideanPoint> tmpArm(armID, groupPoints[i], groupPoints[j]);
            groupIDtoArmID[std::make_pair(i, j)] = armID;
            armsVec.push_back(tmpArm);
            armID++;
        }
    }

    std::vector<std::shared_ptr<SquaredEuclideanPoint> > groupPointTmp1;
    std::vector<std::shared_ptr<SquaredEuclideanPoint> > groupPointTmp2;
    groupPointTmp1.push_back(sharedPtrPointsVec[1]);
    groupPointTmp2.push_back(sharedPtrPointsVec[0]);
    GroupPoint<SquaredEuclideanPoint> gpTmp1(groupPointTmp1, d, 1);
    GroupPoint<SquaredEuclideanPoint> gpTmp2(groupPointTmp2, d, 0);
    ArmHeirarchical<SquaredEuclideanPoint> tmpArm(0, gpTmp1, gpTmp2);
    groupIDtoArmID[std::make_pair(1, 0)] = armID;
    armsVec.push_back(tmpArm);

    UCBDynamic<ArmHeirarchical<SquaredEuclideanPoint> > UCB1(armsVec, delta, 1, 0);

    //MAB initialise
    UCB1.initialise();
    std::cout << "initialised " << std::endl;

    // Step 2: Run clustering

    UCB1.runUCB(10000000);
    std::cout << UCB1.topKArms[0].leftGroupID<< " " << UCB1.topKArms[0].rightGroupID << std::endl;
    for(unsigned long index(0); index < maxPointId ; index++){

    }


//    for (unsigned long i(0); i < 1; i++) {
//
//        //Find the best group points to join
//        ArmHeirarchical<SquaredEuclideanPoint> bestArm = UCB1.bruteBestArm();
////        UCB1.runUCB(1000);
//        std::cout << bestArm.leftGroupID << " " << bestArm.rightGroupID << std::endl;
//    }
}
