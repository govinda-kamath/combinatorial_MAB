//
// Created by Govinda Kamath on 2/12/18.
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

#define  UCB

int main()
{

//    std::string nameConfig = argv[1];
//    long startIndex(atol(argv[2])); // Start index
//    long endIndex(atol(argv[3])); // End index

//     For debugging mode in CLion
    std::string nameConfig = "/Users/govinda/Code/combinatorial_MAB/nominal1.ini";

    // Parameters
    INIReader reader(nameConfig);
    if (reader.ParseError() < 0) {
        std::cout << "Can't load "<< nameConfig << std::endl;
        return 1;
    }
    std::string directoryPath = reader.Get("path", "directory", "");
    std::string saveFilePath =reader.Get("path", "saveFilePath", "test.output");
    std::string fileSuffix = reader.Get("path", "suffix", "");
    unsigned numberOfInitialPulls = (unsigned) reader.GetInteger("UCB", "numberOfInitialPulls_knn", 100);
    float delta = (float) reader.GetReal("UCB", "delta", 0.1);
    std::vector<std::string>  pathsToImages;
    std::vector<SquaredEuclideanPoint> pointsVec;
    std::vector<std::shared_ptr<SquaredEuclideanPoint>> sharedPtrPointsVec;
    std::vector<GroupPoint<SquaredEuclideanPoint>  > groupPoints;
    std::vector<std::string> groupPointsNames;
    std::vector<ArmHeirarchical<SquaredEuclideanPoint> > armsVec;
    // Stores map from group id to arm id. Used to removes "irrelevant" arms
    std::unordered_map<std::pair<unsigned long, unsigned long >, unsigned long, utils::pair_hash> groupIDtoArmID;

    //Maintains the valid group points currently.
    //Used while creating new arms
    std::unordered_set<unsigned  long > groupPointsInPlay;


    utils::getPathToFile(pathsToImages, directoryPath, fileSuffix);
    utils::vectorsToPoints(pointsVec, pathsToImages);
    // Step 1: Initialize all the leaves
    unsigned long armID = 0;
    unsigned long n = pointsVec.size();
    unsigned long d = pointsVec[0].getVecSize();
    unsigned long maxGroupPointId (0);
    std::cout << "Running Hierarchical clustering for " << n << " points" << std::endl;

    groupIDtoArmID.reserve(n*n);
    groupPoints.reserve(n*n);

    for(unsigned long i(0); i< n*n; i++){
        groupPoints.push_back(GroupPoint<SquaredEuclideanPoint> ()); //Todo: Bad Code
        groupPointsNames.push_back("");
    }

    for (unsigned long i(0); i < n; i++) {
        sharedPtrPointsVec.push_back(std::make_shared<SquaredEuclideanPoint>(pointsVec[i]));
        std::vector<std::shared_ptr<SquaredEuclideanPoint> > groupPointTmp1;
        groupPointTmp1.push_back(sharedPtrPointsVec[i]);
        GroupPoint<SquaredEuclideanPoint> gpTmp1(groupPointTmp1, d, maxGroupPointId);
        groupPoints[i] = gpTmp1;
        groupPointsNames[i] =  std::to_string(i);
        groupPointsInPlay.insert(maxGroupPointId);
        maxGroupPointId++;

    }

    for (unsigned long i(0); i < n; i++) {
        for (unsigned long j(0); j < i; j++) {
            ArmHeirarchical<SquaredEuclideanPoint> tmpArm(armID, groupPoints[i], groupPoints[j]);
            groupIDtoArmID[std::make_pair(i,j)] = armID;
            armsVec.push_back(tmpArm);
            armID ++;
        }
    }

    UCBDynamic<ArmHeirarchical<SquaredEuclideanPoint> > UCB1(armsVec, delta, 1, 0);
    UCB1.initialise(100);
    //Brute Only
    UCB1.armsKeepFromArmsContainerBrute();
    std::chrono::system_clock::time_point timeStartStart = std::chrono::system_clock::now();

    // Step 2: Run clustering
    for(unsigned long i(0); i < n-2; i++){
        std::chrono::system_clock::time_point timeStart = std::chrono::system_clock::now();
        //Find the best group points to join
        ArmHeirarchical<SquaredEuclideanPoint> bestArm = UCB1.bruteBestArm();


        std::cout << "Step " << i
                << " Arm ID " << bestArm.id;

        std::cout << " Avg Pulls" << d;

        auto left  = bestArm.leftGroupID;
        auto right = bestArm.rightGroupID;
//        if(i<10){
//            std::cout<< pathsToImages[left] << std::endl;
//            std::cout<< pathsToImages[right] << std::endl;
//        }
        std::cout <<" Left " << groupPointsNames[left]<< " Right " << groupPointsNames[right];

        for(unsigned long index(0); index < maxGroupPointId ; index++){
            if(groupIDtoArmID.find(std::make_pair(left, index)) != groupIDtoArmID.end()){
                UCB1.markForRemoval(groupIDtoArmID[std::make_pair(left, index)]);
            }
            if(groupIDtoArmID.find(std::make_pair(index, left)) != groupIDtoArmID.end()){
                UCB1.markForRemoval(groupIDtoArmID[std::make_pair(index, left)]);
            }
            if(groupIDtoArmID.find(std::make_pair(right, index)) != groupIDtoArmID.end()){
                UCB1.markForRemoval(groupIDtoArmID[std::make_pair(right, index)]);
            }
            if(groupIDtoArmID.find(std::make_pair(index, right)) != groupIDtoArmID.end()){
                UCB1.markForRemoval(groupIDtoArmID[std::make_pair(index, right)]);
            }
        }

        groupPointsInPlay.erase(left);
        groupPointsInPlay.erase(right);


        //Create the vector of pointers to the points in the new group point
        std::vector<std::shared_ptr<SquaredEuclideanPoint > > newGroupPoint;
        newGroupPoint = bestArm.leftGroupPoint->groupPoint;
        newGroupPoint.insert(newGroupPoint.end(), bestArm.rightGroupPoint->groupPoint.begin(),
                                                bestArm.rightGroupPoint->groupPoint.end());


        //create a new group point
        GroupPoint<SquaredEuclideanPoint> gpTmp(newGroupPoint, d, maxGroupPointId);
        groupPoints[maxGroupPointId] = gpTmp;
        groupPointsNames[maxGroupPointId] = groupPointsNames[left]+" "+groupPointsNames[right];

        unsigned  long idOfInsertedPoint = maxGroupPointId;
        maxGroupPointId++;

        //Add and initialise best arm
        for (const auto& pointID: groupPointsInPlay) {
            ArmHeirarchical<SquaredEuclideanPoint> tmpArm(armID, groupPoints[idOfInsertedPoint], groupPoints[pointID]);
            groupIDtoArmID[std::make_pair(idOfInsertedPoint, pointID)] = armID;

            UCB1.initialiseAndAddNewArmBrute(tmpArm); //Brute method
            armID ++;
        }

        //Adding inserted point to points in play
        groupPointsInPlay.insert(idOfInsertedPoint);
        std::chrono::system_clock::time_point timeEnd = std::chrono::system_clock::now();
        long long int trueMeanTime = std::chrono::duration_cast<std::chrono::milliseconds>
                    (timeEnd-timeStart).count();
        std::cout << " Time = " << trueMeanTime << "(ms)" << std::endl;
    }

    std::chrono::system_clock::time_point timeEndEnd = std::chrono::system_clock::now();
    long long int totalTime = std::chrono::duration_cast<std::chrono::milliseconds>
            (timeEndEnd-timeStartStart).count();
    std::cout << "Total Time = " << totalTime << "(ms)" << std::endl;
}