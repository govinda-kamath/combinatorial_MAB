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


int main(int argc, char *argv[])
{
    std::string nameConfig = argv[1];

//     For debugging mode in CLion
//    std::string nameConfig = "/Users/vivekkumarbagaria/Code/combinatorial_MAB/nominal.ini";

    // Parameters
    INIReader reader(nameConfig);
    if (reader.ParseError() < 0) {
        std::cout << "Can't load "<< nameConfig << std::endl;
        return 1;
    }

    // Loading Hyper parameters and data sizes
    std::string directoryPath = reader.Get("path", "h_directory", "");
    std::string saveFilePath = reader.Get("path", "saveFilePathHeirarchical", "test.output");
    std::string fileSuffix = reader.Get("path", "suffix", "");
    unsigned numberOfInitialPulls = (unsigned) reader.GetInteger("UCB", "numberOfInitialPulls_knn", 100);
    unsigned sampleSize = (unsigned) reader.GetInteger("UCB", "sampleSize", 32);
    float delta = (float) reader.GetReal("UCB", "delta", 0.1);


    std::vector<std::string>  pathsToImages;
    std::vector<SquaredEuclideanPoint> pointsVec;
    /*
     * Since vectors are shared across arms, they use shared pointers
     * to share the memory location of vectors. This is pretty obvious but very
     * important
     */
    std::vector<std::shared_ptr<SquaredEuclideanPoint> > sharedPtrPointsVec;

    std::vector<GroupPoint<SquaredEuclideanPoint> > groupPoints;
    std::vector<std::string> groupPointsNames; // Only for debugging purpose: Stores the groups in text format
    std::vector<ArmHeirarchical<SquaredEuclideanPoint> > armsVec;
    // Stores map from group id to arm id. Used to removes "irrelevant" arms
    std::unordered_map<std::pair<unsigned long, unsigned long >, unsigned long, utils::pair_hash> groupIDtoArmID;

    //Maintains the valid group points currently.
    //Used while creating new arms
    std::unordered_set<unsigned long > activeGroups;

    std::ofstream saveFile;
    saveFile.open (saveFilePath, std::ofstream::out | std::ofstream::trunc);
    std::ofstream saveFileBrute;
    saveFileBrute.open (saveFilePath+"brute", std::ofstream::out | std::ofstream::trunc);


    utils::getPathToFile(pathsToImages, directoryPath, fileSuffix); // Loads the filepaths into an array
    utils::vectorsToPoints(pointsVec, pathsToImages); //Loads the images from the locations in above array to pointsVec

    // Step 1: Initialize all the leaves
    unsigned long armID = 0;
    unsigned long n = 200; // pointsVec.size();
    unsigned long d = pointsVec[0].getVecSize();
    unsigned long groupPointId (0);
    std::cout << "Running Hierarchical clustering for " << n << " points" << std::endl;

    groupIDtoArmID.reserve(n*n);
    groupPoints.reserve(n*n);

    // Pushing null objects
    for(unsigned long i(0); i< n*n; i++){
        groupPoints.push_back(GroupPoint<SquaredEuclideanPoint> ()); //Todo: Bad Code
        groupPointsNames.push_back("");
    }

    for (unsigned long i(0); i < n; i++) {
        sharedPtrPointsVec.push_back(std::make_shared<SquaredEuclideanPoint>(pointsVec[i]));
        std::vector<std::shared_ptr<SquaredEuclideanPoint> > groupPointTmp1;
        groupPointTmp1.push_back(sharedPtrPointsVec[i]);
        GroupPoint<SquaredEuclideanPoint> gpTmp1(groupPointTmp1, d, groupPointId);
        groupPoints[i] = gpTmp1;
        groupPointsNames[i] =  std::to_string(i);
        //Updating the gropu gorup to the active list
        activeGroups.insert(groupPointId);
        groupPointId++;
    }

    // Adding the first set of n choose 2 arms.
    for (unsigned long i(0); i < n; i++) {
        for (unsigned long j(0); j < i; j++) {
            ArmHeirarchical<SquaredEuclideanPoint> tmpArm(armID, groupPoints[i], groupPoints[j]);
//            float lore  = tmpArm.trueMean();
            groupIDtoArmID[std::make_pair(i,j)] = armID;
            armsVec.push_back(tmpArm);
            armID++;
        }
    }

    UCBDynamic<ArmHeirarchical<SquaredEuclideanPoint> > UCB1(armsVec, delta, 1, 0, sampleSize);

#define UCB
#ifdef UCB
    std::cout << "Running MAB with sample size " << sampleSize << std::endl;
    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
    UCB1.initialise(numberOfInitialPulls);
    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
    long long int Tt = std::chrono::duration_cast<std::chrono::milliseconds>
            (t2-t1).count();
    std::cout << "Init Time = " << Tt << "(ms)" << std::endl;

#else
    std::cout << "Running Brute" << std::endl;
    UCB1.armsKeepFromArmsContainerBrute();
#endif

    std::chrono::system_clock::time_point timeStartStart = std::chrono::system_clock::now();
    unsigned long prevGlobalNumberPulls = UCB1.globalNumberOfPulls;

    // Step 2: Run clustering
    for(unsigned long i(0); i < n-2; i++){
        std::chrono::system_clock::time_point timeStart = std::chrono::system_clock::now();
        //Find the best group points to join
#ifdef UCB
        UCB1.runUCB(n*d);
        ArmHeirarchical<SquaredEuclideanPoint> bestArm = UCB1.topKArms.back();
#else
        ArmHeirarchical<SquaredEuclideanPoint> bestArm = UCB1.bruteBestArms()[0];
#endif

//        std::cout << "Step " << i  << " Arm ID " << bestArm.id  << " Avg Pulls" << UCB1.globalNumberOfPulls/(n*n*.5);
        auto left  = bestArm.leftGroupID;
        auto right = bestArm.rightGroupID;


//        std::cout <<" Left " << groupPointsNames[left]<< " Right " << groupPointsNames[right];

        //Removing arms containing either of the above left and right group of points.
        for(unsigned long index(0); index < groupPointId ; index++){
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

        // Removing left and right groups from the list of groups (nodes)
        activeGroups.erase(left);
        activeGroups.erase(right);

        //Creating a new vector of pointers by combining left and right group points
        std::vector<std::shared_ptr<SquaredEuclideanPoint > > newGroupPoint;
        newGroupPoint = bestArm.leftGroupPoint->groupPoint;
        newGroupPoint.insert(newGroupPoint.end(), bestArm.rightGroupPoint->groupPoint.begin(),
                             bestArm.rightGroupPoint->groupPoint.end());


        //Create a new group point from the new vector of pointers
        GroupPoint<SquaredEuclideanPoint> newCombinedGroupPoint(newGroupPoint, d, groupPointId);
        groupPoints[groupPointId] = newCombinedGroupPoint;
        groupPointsNames[groupPointId] = groupPointsNames[left]+" "+groupPointsNames[right]; //Debug

        /*
         * Remove old arms, add and initialise new arms.
         * For each active group, add an arm comprising of that group and the new group created above
         */
        for (const auto& pointID: activeGroups) {
            ArmHeirarchical<SquaredEuclideanPoint> tmpArm(armID, groupPoints[groupPointId], groupPoints[pointID]);
            groupIDtoArmID[std::make_pair(groupPointId, pointID)] = armID;

            unsigned long leftArmRemovedId;
            unsigned long rightArmRemovedId;

            // Removing arms
            if(groupIDtoArmID.find(std::make_pair(left, pointID)) != groupIDtoArmID.end()){
                leftArmRemovedId = groupIDtoArmID[std::make_pair(left, pointID)];
            }
            else if (groupIDtoArmID.find(std::make_pair(pointID, left)) != groupIDtoArmID.end()){
                leftArmRemovedId = groupIDtoArmID[std::make_pair(pointID, left)];
            }
            else{
                throw std::runtime_error("[Unexpected behaviour]: Marked left arm's id not found.");
            }

            if(groupIDtoArmID.find(std::make_pair(right, pointID)) != groupIDtoArmID.end()){
                rightArmRemovedId = groupIDtoArmID[std::make_pair(right, pointID)];
            }
            else if (groupIDtoArmID.find(std::make_pair(pointID, right)) != groupIDtoArmID.end()){
                rightArmRemovedId = groupIDtoArmID[std::make_pair(pointID, right)];
            }
            else{
                throw std::runtime_error("[Unexpected behaviour]: Marked right arm's id not found.");
            }

            //Combining the data from the two arms and shift that data to the new arm
            unsigned long numArmPulls(0);
            float armSumOfPulls(0.0);
            float armSumOfSquaresOfPulls(0.0);
            float trueMeanValue(0.0);

            long minNumberOfPulls = std::min(UCB1.armStates[leftArmRemovedId].numberOfPulls,
                                             UCB1.armStates[rightArmRemovedId].numberOfPulls);



            float leftSize  = bestArm.leftGroupPoint->noOfPoints;
            float rightSize = bestArm.rightGroupPoint->noOfPoints;

            float f = (leftSize)/(leftSize+rightSize);
            unsigned long n1 = UCB1.armStates[leftArmRemovedId].numberOfPulls;
            unsigned long n2 = UCB1.armStates[rightArmRemovedId].numberOfPulls;
            float mu1 = UCB1.armStates[leftArmRemovedId].sumOfPulls/n1;
            float mu2 = UCB1.armStates[rightArmRemovedId].sumOfPulls/n2;

            float s1 = UCB1.armStates[leftArmRemovedId].sumOfSquaresOfPulls/n1;
            float s2 = UCB1.armStates[rightArmRemovedId].sumOfSquaresOfPulls/n2;
            float sigma1 = std::sqrt(s1-mu1*mu1);
            float sigma2 = std::sqrt(s2-mu2*mu2);

            float t1 = UCB1.armStates[leftArmRemovedId].trueMeanValue;
            float t2 = UCB1.armStates[rightArmRemovedId].trueMeanValue;

            float uEff = f*mu1+(1-f)*mu2;
            float sEff = f*s1+(1-f)*s2;
            float sigmaEff = std::sqrt(sEff - uEff*uEff);

            numArmPulls = std::min(n1,n2);
            armSumOfPulls = uEff*numArmPulls;
            armSumOfSquaresOfPulls = sEff*numArmPulls;
            trueMeanValue = f*t1+(1-f)*t2;

//            std::cout << "IMMMP " << groupPointsNames[pointID];
//            std::cout << "\t" << groupPointsNames[groupPointId];
//            std::cout << "\t" << groupPointsNames[left];
//            std::cout << "\t" << groupPointsNames[right];
//            std::cout << "\t" << t1;
//            std::cout << "\t" << t2;
//            std::cout << "\t" << trueMeanValue << std::endl;
//

//            std::cout << "IMMMP " << groupPointsNames[pointID];
//            std::cout << "\t" << groupPointsNames[groupPointId];
//            std::cout << "\t" << n1 << " " << n2;
//            std::cout << "\t" << mu1 << " " << t1;
//            std::cout << "\t" << mu2 << " " << t2;
//            std::cout << "\t" << uEff << " " << trueMeanValue << std::endl;
////
//            if ( (mu1-sigma1 > t1) or(mu1+sigma1 < t1)){
//                std::cout << "Point1 group!! ";
//                std::cout << "\t" << leftArmRemovedId;
//                std::cout << "Left: " << mu1-sigma1 << " " << " " << mu1 <<"," << t1 << " " <<mu1+sigma1;
//
//                bad++;
//
//            }
//            if ( (mu2-sigma2 > t2) or(mu2+sigma2 < t2)){
//                std::cout << "\nPoint2 Wrong!! ";
//                std::cout << "\t" << rightArmRemovedId;
//                std::cout << "\nRight: " << mu2-sigma2 << " " << mu2 <<"," << t2 << " " <<mu2+sigma2;
//
//            std::cout << "IMMMP " << groupPointsNames[pointID];
//            std::cout << "\t" << groupPointsNames[groupPointId];
//            std::cout << "\t" << groupPointsNames[left];
//            std::cout << "\t" << groupPointsNames[right];
//            std::cout << "\t" << t1;
//            std::cout << "\t" << t2;
//            std::cout << "\t" << trueMeanValue << std::endl;
//            bad++;
//            }
//            if ( (uEff-sigmaEff> trueMeanValue) or(uEff+sigmaEff< trueMeanValue) ){
//                std::cout << "\nPoint Combination Wrong!!";
//                std::cout << "\nCombination: " << uEff-sigmaEff << " " << uEff <<"," << trueMeanValue << " " <<uEff+sigmaEff;
//            }
//            std::cout << "Left: " << mu1-sigma1 << " " << " " << mu1 <<"," << t1 << " " <<mu1+sigma1;
//            std::cout << "\nRight: " << mu2-sigma2 << " " << mu2 <<"," << t2 << " " <<mu2+sigma2;
//            std::cout << "\nCombination: " << u << mu1 <<"," << t1 << " " <<mu1+sigma1;
//            std::cout << "\nRight: " << mu2-sigma2 << " " << mu2 <<"," << t2 << " " <<mu2+sigma2;
//            std::cout << "\nCombination: " << uEff-sigmaEff << " " << uEff <<"," << trueMeanValue << " " <<uEff+sigmaEff;

//            std::cout << "ratio " << f <<" n1 " << n1 << " n2 " << n2;
////            std::cout << " mu1 " << mu1 << " mu2 " <<mu2;
////            std::cout << " s1 " << s1 << " s2 " <<s2;
////            std::cout <<" Left Size " << leftSize << " Right Size " <<rightSize;
//            std::cout << "LG " << uEff - sigmaEff;
//            std::cout <<" uEff " << uEff;
//            std::cout << "UG " << uEff + sigmaEff;
//            std::cout <<" TrueMean " << trueMeanValue;
////            std::cout << "\t" << groupPointsNames[leftArmRemovedId];
////            std::cout << "\t" << groupPointsNames[rightArmRemovedId];
//            std::cout << std::endl;

//            if(UCB1.armStates.find(leftArmRemovedId) != UCB1.armStates.end()){
//                numArmPulls += UCB1.armStates[leftArmRemovedId].numberOfPulls;
//                armSumOfPulls += UCB1.armStates[leftArmRemovedId].sumOfPulls;
//                armSumOfSquaresOfPulls += UCB1.armStates[leftArmRemovedId].sumOfSquaresOfPulls;
//                trueMeanValue += UCB1.armStates[leftArmRemovedId].trueMeanValue;
//            }
//            else{
//                throw std::runtime_error("[Unexpected behaviour]: Marked left arm's state not found.");
//            }
//            if(UCB1.armStates.find(rightArmRemovedId) != UCB1.armStates.end()){
//                numArmPulls += UCB1.armStates[rightArmRemovedId].numberOfPulls;
//                armSumOfPulls += UCB1.armStates[rightArmRemovedId].sumOfPulls;
//                armSumOfSquaresOfPulls += UCB1.armStates[rightArmRemovedId].sumOfSquaresOfPulls;
//                trueMeanValue += UCB1.armStates[rightArmRemovedId].trueMeanValue;
//            }
//            else{
//                throw std::runtime_error("[Unexpected behaviour]: Marked right arm's state not found.");
//            }

//            std::cout << "Left pulls " << UCB1.armStates[leftArmRemovedId].numberOfPulls;
//            std::cout << " Right pulls " << UCB1.armStates[rightArmRemovedId].numberOfPulls;
//            std::cout << std::endl;

            //Create a new arm
#ifdef UCB
//            tmpArm.warmInitialise(numArmPulls, armSumOfPulls, armSumOfSquaresOfPulls, trueMeanValue);
//            UCB1.initialiseAndAddNewArm(tmpArm, 0); //0 initial pulls
            tmpArm.warmInitialise(0, 0, 0, INFINITY);
            UCB1.initialiseAndAddNewArm(tmpArm, numberOfInitialPulls); //0 initial pulls
#else
            tmpArm.trueMeanValue = trueMeanValue;
            tmpArm.lowerConfidenceBound = tmpArm.trueMeanValue;
            tmpArm.upperConfidenceBound = tmpArm.trueMeanValue;
            tmpArm.estimateOfMean = tmpArm.trueMeanValue;
            UCB1.initialiseAndAddNewArmBrute(tmpArm, 0);
#endif
            armID++;
        }

        //Adding inserted group point in the active group
        activeGroups.insert(groupPointId);
        groupPointId++;

        std::chrono::system_clock::time_point timeEnd = std::chrono::system_clock::now();
        long long int trueMeanTime = std::chrono::duration_cast<std::chrono::milliseconds>
                (timeEnd-timeStart).count();

        std::chrono::system_clock::time_point timeEndEnd = std::chrono::system_clock::now();
        long long int totalTime = std::chrono::duration_cast<std::chrono::milliseconds>
                (timeEndEnd-timeStartStart).count();

        float localVar = bestArm.estimateOfSecondMoment - std::pow(bestArm.estimateOfMean,2);
        std::cout   << "Step " << i
                << "\tTrue. = " << bestArm.trueMean()
//                << "\tSecond True. = " << UCB1.topValidArm().trueMean()
                << "\tEst. = " << bestArm.estimateOfMean
                << "\tId. 1= " << bestArm.leftGroupID
                << "\tId. 2= " << bestArm.rightGroupID
//                << "\tSecond Id. 1= " << UCB1.topValidArm().leftGroupID
//                << "\tSecond Id. 2= " << UCB1.topValidArm().rightGroupID
//                    << "\tTime = " << trueMeanTime
//                    << "\tT-Time = " << totalTime
                  << "\tGlobal number of Pulls = " << UCB1.globalNumberOfPulls
//                    << "\tNew no. of Puls = " << UCB1.globalNumberOfPulls - prevGlobalNumberPulls
////                  << "\ttime/Pullslaststep = " << (float)trueMeanTime/(UCB1.globalNumberOfPulls - prevGlobalNumberPulls)
//                    << "\ttime/ Pullstotal = " << (float)totalTime/(UCB1.globalNumberOfPulls)
//                    << "\tGlobal Sig= " << UCB1.globalSigma
//                    << "\tActive group size=" << activeGroups.size()
//                << "\tBest Arm local Sigma= " << localVar
//                << "\tBest Arm local estimate= " << bestArm.estimateOfMean
//                << "\tBest Arm local var= " << bestArm.estimateOfSecondMoment
//                << "\tLeft= " << groupPointsNames[left]
//                    << "\tRight= " << groupPointsNames[right]
                    << std::endl;
        prevGlobalNumberPulls = UCB1.globalNumberOfPulls;
        saveFile <<  groupPointId-1 << "," <<  left << "," << right << std::endl;

#ifndef UCB
        saveFileBrute <<  groupPointId-1 << "," <<  left << "," << right << std::endl;
#endif
    }
    saveFile <<  groupPointId;
    std::cout << "GroupID\t" << groupPointId << std::endl;
    for (const auto& pointID: activeGroups) {
        saveFile <<  "," << pointID ;
        std::cout << "Active left\t" << pointID << std::endl;
    }
    saveFile << std::endl;
    std::chrono::system_clock::time_point timeEndEnd = std::chrono::system_clock::now();
    long long int totalTime = std::chrono::duration_cast<std::chrono::milliseconds>
            (timeEndEnd-timeStartStart).count();
    std::cout << "Total Time = " << totalTime << "(ms)"
            << "Global number of Pulls = " << UCB1.globalNumberOfPulls/(0.5*n*n)
            << "Global Sigma= " << UCB1.globalSigma
            << std::endl;

#ifdef UCB
    std::cout<< "Average distances evaluated " << UCB1.globalNumberOfPulls/(0.5*n*n) << std::endl;
#else
    std::cout<< "Average distances evaluated " << d << std::endl;
#endif
}
