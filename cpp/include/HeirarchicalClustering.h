//
// Created by Vivek Kumar Bagaria on 2/3/18.
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
#include <dlib/image_io.h>
#include <dlib/image_transforms.h>
#include "Points.h"
#include "Arms.h"
#include "UCB_dynamic.h"
#include <stdexcept>
#include "utils.h"

/*
 * Finds the nearest k neighbhours of pointsVectorLeft in pointsVectorRight
 */
template <class templatePoint>
class Heirarchical{
public:

    std::vector<templatePoint> pointsVector;
    unsigned numberOfInitialPulls; // UCB Parameters
    float delta; // UCB Parameters
    unsigned int sampleSize; //UCB Parameters
    unsigned long n, d, armID, groupPointId;
    char algo; //mab or brute?


    //Variables used
    /*
     * Since vectors are shared across arms, we use shared pointers to share
     * the memory location of vectors. This is pretty obvious but very important
     */
    std::vector<std::shared_ptr<templatePoint> > sharedPtrPointsVec;

    std::vector<GroupPoint<templatePoint> > groupPoints;
    std::vector<std::string> groupPointsNames; // Only for debugging purpose: Stores the groups in text format
    std::vector<ArmHeirarchical<templatePoint> > armsVec;
    // Stores map from group id to arm id. Used to removes "irrelevant" arms
    std::unordered_map<std::pair<unsigned long, unsigned long>, unsigned long, utils::pair_hash> groupIDtoArmID;
    //Maintains the valid group points currently and used while creating new arms
    std::unordered_set<unsigned long> activeGroups;


    //Outputs Variables
    std::vector<std::vector<unsigned int>> merges;
    std::vector<ArmKNN<templatePoint>> mergesBrute;
    std::unordered_map<unsigned long, unsigned long> finalNumberOfPulls;
    std::vector<unsigned long> finalSortedOrder;
    float initTime;
    float runTime;
    std::ofstream saveFile; //Stats output file
    std::ofstream graphSaveFile; //Graph output file

    Heirarchical( const std::vector<templatePoint> &pVec, unsigned noOfInitialPulls, float deltaAccuracy,
        unsigned int sSize, std::string sFilePath, std::string gFilePath, char alg, unsigned long nn) {

        pointsVector = pVec;
        d = pointsVector[0].getVecSize();
        numberOfInitialPulls = noOfInitialPulls;
        delta = deltaAccuracy;
        sampleSize = sSize;
        n = nn;
        algo = alg;
        gFilePath = gFilePath+"n_"+std::to_string(n)+"_d_"+std::to_string(d);
        sFilePath = sFilePath+"n_"+std::to_string(n)+"_d_"+std::to_string(d);

        if (algo == 'm') {
            graphSaveFile.open(gFilePath, std::ofstream::out | std::ofstream::trunc);
            saveFile.open(sFilePath, std::ofstream::out | std::ofstream::trunc);
        } else {
            graphSaveFile.open(gFilePath + "brute", std::ofstream::out | std::ofstream::trunc);
            saveFile.open(sFilePath + "brute", std::ofstream::out | std::ofstream::trunc);
        }

        setup(); //Todo: Check
    }



    void setup(){

        armID = 0;
        groupPointId = 0;
        groupIDtoArmID.reserve(n * n);
        groupPoints.reserve(n * n);

        // Pushing null objects
        std::cout << "Creating n^2 Group points";
        for (unsigned long i(0); i < n * n; i++) {
            groupPoints.push_back(GroupPoint<templatePoint>()); //Todo: Bad Code
            groupPointsNames.push_back("");
        }

        for (unsigned long i(0); i < n; i++) {
            sharedPtrPointsVec.push_back(std::make_shared<templatePoint>(pointsVector[i]));
            std::vector<std::shared_ptr<templatePoint> > groupPointTmp1;
            groupPointTmp1.push_back(sharedPtrPointsVec[i]);
            GroupPoint<templatePoint> gpTmp1(groupPointTmp1, d, groupPointId);
            groupPoints[i] = gpTmp1;
            groupPointsNames[i] = std::to_string(i);
            //Updating the  group to the active list
            activeGroups.insert(groupPointId);
            groupPointId++;
        }

        std::cout << "Adding the first n*n/2 arms" << std::endl;
        // Adding the first set of n choose 2 arms.
        for (unsigned long i(0); i < n; i++) {
            if (i % (int) std::pow(n, .5) == 0) {
                std::cout << i << " arms added" << std::endl;
            }
            for (unsigned long j(0); j < i; j++) {
                ArmHeirarchical<templatePoint> tmpArm(armID, groupPoints[i], groupPoints[j]);
                groupIDtoArmID[std::make_pair(i, j)] = armID;
                armsVec.push_back(tmpArm);
                armID++;
            }
        }

    }

    void run(){

        UCBDynamic<ArmHeirarchical<templatePoint> > UCB1(armsVec, delta, 1, 0, sampleSize);

        if (algo == 'm') {
            std::cout << "Running MAB" << std::endl;
            UCB1.initialise(numberOfInitialPulls);

        } else {
            std::cout << "Running Brute" << std::endl;
            UCB1.armsKeepFromArmsContainerBrute();
        }
        std::cout << "Init done "<< std::endl;

        std::chrono::system_clock::time_point timeStartStart = std::chrono::system_clock::now();

        // Step 2: Run clustering
        for (unsigned long i(0); i < n - 2; i++) {

            //Find the best group points to join
            ArmHeirarchical<templatePoint> *bestArm;
            if (algo == 'm') {
                UCB1.runUCB(n * d);
                bestArm = &UCB1.topKArms.back();
            } else {
                std::vector<ArmHeirarchical<templatePoint>> bestArms = UCB1.bruteBestArms();
                bestArm = &bestArms[0];
            }
            auto left = bestArm->leftGroupID;
            auto right = bestArm->rightGroupID;

            // Removing left and right groups from the list of groups (nodes)
            activeGroups.erase(left);
            activeGroups.erase(right);
            //Removing arms containing either of the above left and right group of points.
            for (unsigned long index(0); index < groupPointId; index++) {

                if (groupIDtoArmID.find(std::make_pair(left, index)) != groupIDtoArmID.end()) {
                    UCB1.markForRemoval(groupIDtoArmID[std::make_pair(left, index)]);
                }
                if (groupIDtoArmID.find(std::make_pair(index, left)) != groupIDtoArmID.end()) {
                    UCB1.markForRemoval(groupIDtoArmID[std::make_pair(index, left)]);
                }
                if (groupIDtoArmID.find(std::make_pair(right, index)) != groupIDtoArmID.end()) {
                    UCB1.markForRemoval(groupIDtoArmID[std::make_pair(right, index)]);
                }
                if (groupIDtoArmID.find(std::make_pair(index, right)) != groupIDtoArmID.end()) {
                    UCB1.markForRemoval(groupIDtoArmID[std::make_pair(index, right)]);
                }
            }


            //Creating a new group by combining left and right group points
            std::vector<std::shared_ptr<templatePoint> > newGroupPoint;
            newGroupPoint = bestArm->leftGroupPoint->groupPoint;
            newGroupPoint.insert(newGroupPoint.end(), bestArm->rightGroupPoint->groupPoint.begin(),
                                 bestArm->rightGroupPoint->groupPoint.end());
            GroupPoint<templatePoint> newCombinedGroupPoint(newGroupPoint, d, groupPointId);
            groupPoints[groupPointId] = newCombinedGroupPoint;
            groupPointsNames[groupPointId] = groupPointsNames[left] + " " + groupPointsNames[right]; //Debug


            /*
             * Remove old arms, add and initialise new arms.
             * For each active group, add an arm comprising of that group and the new group created above
             */
            for (const auto &pointID: activeGroups) {

                ArmHeirarchical<templatePoint> tmpArm(armID, groupPoints[groupPointId], groupPoints[pointID]);
                groupIDtoArmID[std::make_pair(groupPointId, pointID)] = armID;

                unsigned long leftArmRemovedId;
                unsigned long rightArmRemovedId;

                // Removing arms
                if (groupIDtoArmID.find(std::make_pair(left, pointID)) != groupIDtoArmID.end()) {
                    leftArmRemovedId = groupIDtoArmID[std::make_pair(left, pointID)];
                } else if (groupIDtoArmID.find(std::make_pair(pointID, left)) != groupIDtoArmID.end()) {
                    leftArmRemovedId = groupIDtoArmID[std::make_pair(pointID, left)];
                } else {
                    throw std::runtime_error("[Unexpected behaviour]: Marked left arm's id not found.");
                }

                if (groupIDtoArmID.find(std::make_pair(right, pointID)) != groupIDtoArmID.end()) {
                    rightArmRemovedId = groupIDtoArmID[std::make_pair(right, pointID)];
                } else if (groupIDtoArmID.find(std::make_pair(pointID, right)) != groupIDtoArmID.end()) {
                    rightArmRemovedId = groupIDtoArmID[std::make_pair(pointID, right)];
                } else {
                    throw std::runtime_error("[Unexpected behaviour]: Marked right arm's id not found.");
                }


                //Create a new arm
                if (algo == 'm') {
                    tmpArm.warmInitialise(0, 0, 0, INFINITY);
                    UCB1.initialiseAndAddNewArm(tmpArm, numberOfInitialPulls);
                } else {
                    float leftSize = bestArm->leftGroupPoint->noOfPoints;
                    float rightSize = bestArm->rightGroupPoint->noOfPoints;
                    float f = (leftSize) / (leftSize + rightSize);
                    float t1 = UCB1.armStates[leftArmRemovedId].trueMeanValue;
                    float t2 = UCB1.armStates[rightArmRemovedId].trueMeanValue;
                    float trueMeanValue =  f * t1 + (1 - f) * t2;
                    tmpArm.trueMeanValue = trueMeanValue;
                    tmpArm.lowerConfidenceBound = tmpArm.trueMeanValue;
                    tmpArm.upperConfidenceBound = tmpArm.trueMeanValue;
                    tmpArm.estimateOfMean = tmpArm.trueMeanValue;
                    UCB1.initialiseAndAddNewArmBrute(tmpArm, 0);
                }
                armID++;
            }

            //Adding inserted group point in the active group
            activeGroups.insert(groupPointId);
            groupPointId++;


            std::cout << "Step " << i
                      << "\tTrue. = " << bestArm->trueMean()
                      //                << "\tSecond True. = " << UCB1.topValidArm().trueMean()
                      << "\tEst. = " << bestArm->estimateOfMean
                      << "\tId. 1= " << bestArm->leftGroupID
                      << "\tId. 2= " << bestArm->rightGroupID
                      << "\tAverage number of Pulls = " << UCB1.globalNumberOfPulls / ((n-1) * (n-1))
                      << std::endl;
            saveFile << "Step " << i
                         << "\tTrue. = " << bestArm->trueMean()
                         << "\tEst. = " << bestArm->estimateOfMean
                         << "\tId. 1= " << bestArm->leftGroupID
                         << "\tId. 2= " << bestArm->rightGroupID
                         << "\tAverage number of Pulls = " << UCB1.globalNumberOfPulls / ((n-1) * (n-1))
                         << std::endl;
            graphSaveFile << groupPointId - 1 << "," << left << "," << right << std::endl;

        }
        graphSaveFile << groupPointId;
        for (const auto &pointID: activeGroups) {
            graphSaveFile << "," << pointID;
//            std::cout << "Active left\t" << pointID << std::endl;
        }
        graphSaveFile << std::endl;
        std::chrono::system_clock::time_point timeEndEnd = std::chrono::system_clock::now();
        long long int totalTime = std::chrono::duration_cast<std::chrono::milliseconds>
                (timeEndEnd - timeStartStart).count();
        std::cout << "Total Time = " << totalTime << "(ms)"
                  << "\nGlobal number of Pulls = " << UCB1.globalNumberOfPulls
                  << "\nGlobal Sigma= " << UCB1.globalSigma
                  << std::endl;
        saveFile   << "TotalTime," << totalTime << std::endl;
        saveFile   << "GlobalnumberofPulls," << UCB1.globalNumberOfPulls << std::endl;
        saveFile   << "GlobalSigma," << UCB1.globalSigma  << std::endl;

        if (algo == 'm') {
            std::cout << "Average distances evaluated " << UCB1.globalNumberOfPulls / ((n-1) * (n-1))
                         << ". Gain = " << d / (UCB1.globalNumberOfPulls / ((n-1) * (n-1))) << std::endl;
            saveFile<< "AveragePulls," << UCB1.globalNumberOfPulls / ((n-1) * (n-1)) << std::endl;

        } else {
            std::cout << "Average distances evaluated " << d << std::endl;
            saveFile<< "AveragePulls," << d << std::endl;
        }
        UCB1.storeExtraTopArms();
        finalNumberOfPulls = UCB1.finalNumberOfPulls;


    }

   void saveAnswers(){
       saveFile << "AllPullsNumber";
       for (unsigned i = 0; i < finalNumberOfPulls.size(); i++) {
           saveFile <<  "," << finalNumberOfPulls[i];
       }
       saveFile << std::endl;
 }
};