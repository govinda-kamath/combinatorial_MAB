#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <random>
#include <thread>
#include <numeric>
#include <map>
#include <dlib/image_io.h>
#include <dlib/image_transforms.h>
#include "Points.h"
#include "INIReader.h"
#include "kmeans.h"
#include "tenXReader.h"


int main(int argc, char *argv[]){

//    std::string nameConfig = argv[1];
    std::string nameConfig = "/Users/vivekkumarbagaria/Code/combinatorial_MAB/nominal.ini";
//
////     Parameters
    INIReader reader(nameConfig);
    if (reader.ParseError() < 0) {
        std::cout << "Can't load " << nameConfig << std::endl;
        return 1;
    }

    // Loading Hyper parameters and data sizes
    unsigned sampleSize = (unsigned) reader.GetInteger("UCB", "sampleSize", 32);
    float delta = (float) reader.GetReal("UCB", "delta", 0.1);
    unsigned numberOfInitialPulls = (unsigned) reader.GetInteger("UCB", "numberOfInitialPulls_knn", 100);
    unsigned long m = (unsigned) reader.GetInteger("UCB", "m", 100);
    unsigned long n = (unsigned) reader.GetInteger("UCB", "n", 100);
    std::string fileName = reader.Get("path", "h5path", "../test_dataset/1M_neurons_matrix_subsampled_20k_5kgenes.h5");

//    delta = delta / (m);

    std::random_device rd;
    std::mt19937 g(9);
    std::normal_distribution<> dist(0, 0.0005);
    std::cout << "Started!" << std::endl;
    std::string str;
    float tmp;

    std::unordered_map<std::pair<unsigned long, unsigned long>, unsigned long, utils::pair_hash> groupIDtoArmID;
    std::vector<std::pair<unsigned long, unsigned long>> armIDtoGroupID;
    //Maintains the valid group points currently and used while creating new arms
    std::unordered_set<unsigned long> verticesAdded;
    std::vector<int> groupMembership;
    unsigned membershipID(0);

//    long n = 10000;
    unsigned long armID = 0;
    std::vector<SquaredEuclideanPoint> dataMatrix;
    std::vector<Arm2DMutualInformation<SquaredEuclideanPoint> > armsVec;

    // Loading 10x data shape
    std::vector<unsigned> shapeData =  tenXReader::get10xMatrixSize(fileName);

    // Loading Sparse matrix
    std::cout << "Reading normalised data" << std::endl;
//    std::vector<std::unordered_map<unsigned long, float> > sparseNormalisedDataMatrix(shapeData[1] );
    unsigned cells(20000), genes(5000);
    std::vector<std::vector<float> > DataMatrix(cells);
    for(unsigned i(0); i< cells ; i ++)
        DataMatrix[i] = std::vector<float>(genes);

    tenXReader::get10xMatrix(fileName, DataMatrix);

    //Arms
    std::vector<SquaredEuclideanPoint> allPointsVec;
    std::cout << "Reading normalised data sparsely part 2" << std::endl;
    utils::vectorsToPoints(allPointsVec, DataMatrix,  n, m);

    armIDtoGroupID.reserve((int)0.5*m*(m-1));
    for(unsigned i(1); i< m ; i ++) {
        for (unsigned j(0); j < i; j++) {
            std::vector<unsigned long> indices = {i, j};
            groupIDtoArmID[std::make_pair(i, j)] = armID;
            armIDtoGroupID.push_back(std::make_pair(i, j));
            std::vector<unsigned long> shuffledRows_(allPointsVec.size());
            std::iota(shuffledRows_.begin(), shuffledRows_.end(), 0);
            std::shuffle(shuffledRows_.begin(), shuffledRows_.end(), g);
            std::chrono::system_clock::time_point sTime = std::chrono::system_clock::now();
            Arm2DMutualInformation<SquaredEuclideanPoint> arm(armID, allPointsVec, indices, shuffledRows_);
            std::chrono::system_clock::time_point eTime = std::chrono::system_clock::now();
            long long int totalTime = std::chrono::duration_cast<std::chrono::milliseconds>
                    (eTime - sTime).count();
            armsVec.push_back(arm);
            armID++;
            std::cout << totalTime << " " << armID << " " << armIDtoGroupID.size()
                      <<  " " << armIDtoGroupID[armID].first
                      << " " << armIDtoGroupID[armID].second  << std::endl;
        }
        groupMembership.push_back(-1);
    }



    std::cout << "UCB go! for steps " << m-1 << "with arms " << armID << std::endl;
    UCBDynamic<Arm2DMutualInformation<SquaredEuclideanPoint> > UCB1(armsVec, delta, 1, 0, sampleSize);
    UCB1.initialise(numberOfInitialPulls);
    for (unsigned long i(0); i < m - 1; i++) {

        //Find the best group points to join
//            ArmHeirarchical<templatePoint> *bestArm;
#ifndef Brute
        UCB1.runUCB(n * 10);
        Arm2DMutualInformation<SquaredEuclideanPoint> bestArm = UCB1.topKArms.back();
        Arm2DMutualInformation<SquaredEuclideanPoint> secondBestArm = UCB1.arms.top();
#else
        float leftSize, rightSize, f;
            std::vector<ArmHeirarchical<templatePoint>> bestArms = UCB1.bruteBestArms();
            ArmHeirarchical<templatePoint> bestArm = bestArms[0];
            leftSize = bestArm.leftGroupPoint->noOfPoints;
            rightSize = bestArm.rightGroupPoint->noOfPoints;
            f = (leftSize) / (leftSize + rightSize);
#endif
        auto left = armIDtoGroupID[bestArm.id].first;
        auto right = armIDtoGroupID[bestArm.id].second;

        auto left2 = armIDtoGroupID[secondBestArm.id].first;
        auto right2 = armIDtoGroupID[secondBestArm.id].second;

        std::cout << "Step " << i
                  //                      << "\tTrue. = " << bestArm.trueMean()
                  //                << "\tSecond True. = " << UCB1.topValidArm().trueMean()
                << "\tLCB " << bestArm.lowerConfidenceBound
                << " Est " << bestArm.estimateOfMean
                << " UCB " << bestArm.upperConfidenceBound
                  << "\tId. 1= " << left
                  << "\tId. 2= " << right
                  << "\tAverage number of Pulls = " << (float) 2*UCB1.globalNumberOfPulls / ((m- 1) * (m - 1))
                << "Pulls " << bestArm.numberOfPulls
                << std::endl;

        std::cout << "Step " << i
                  //                      << "\tTrue. = " << bestArm.trueMean()
                  //                << "\tSecond True. = " << UCB1.topValidArm().trueMean()
                    << "\tLCB " <<secondBestArm.lowerConfidenceBound
                    << " Est " << secondBestArm.estimateOfMean
                    << " UCB " << secondBestArm.upperConfidenceBound
                  << "\tId. 1= " << left2
                  << "\tId. 2= " << right2
                << "\tAverage number of Pulls = " << (float) 2*UCB1.globalNumberOfPulls / ((m- 1) * (m - 1))
                << "Pulls " << secondBestArm.numberOfPulls
                  << std::endl;


        //Removing arms containing either of the above left and right group of points.
        std::cout << "Removing ";



        // Two new vertices
        if(groupMembership[left] == -1 and groupMembership[right] == -1){
            groupMembership[left] = membershipID;
            groupMembership[right] = membershipID;
            membershipID ++;
        }

        else if(groupMembership[left] != -1 and groupMembership[right] == -1){
            for(unsigned index(0); index<m ; index++){
                if (groupMembership[left] == groupMembership[index] and index!=right and index!=left){
                    if (groupIDtoArmID.find(std::make_pair(right, index)) != groupIDtoArmID.end())
                        UCB1.markForRemoval(groupIDtoArmID[std::make_pair(right, index)]);
                    else if (groupIDtoArmID.find(std::make_pair(index, right)) != groupIDtoArmID.end())
                        UCB1.markForRemoval(groupIDtoArmID[std::make_pair(index, right)]);
                    else
                        std::cout << "Problem when removing right, index " << right << " " << index << std::endl;
                }
            }
            groupMembership[right] = groupMembership[left];
        }
        else if(groupMembership[left] == -1 and groupMembership[right] != -1){
            for(unsigned index(0); index<m ; index++){
                if (groupMembership[right] == groupMembership[index] and index!=right and index!=left){
                    if (groupIDtoArmID.find(std::make_pair(left, index)) != groupIDtoArmID.end())
                        UCB1.markForRemoval(groupIDtoArmID[std::make_pair(left, index)]);
                    else if (groupIDtoArmID.find(std::make_pair(index, left)) != groupIDtoArmID.end())
                        UCB1.markForRemoval(groupIDtoArmID[std::make_pair(index, left)]);
                    else
                        std::cout << "Problem when removing right, index " << left << " " << index << std::endl;
                }
            }
            groupMembership[right] = groupMembership[left];
        }

        else {
            for (unsigned index1(0); index1 < m; index1++) {
                for (unsigned index2(0); index2 < m; index2++) {
                    if (groupMembership[index2] == groupMembership[right] and
                        groupMembership[index1] == groupMembership[left]) {
                        if (groupIDtoArmID.find(std::make_pair(index1, index2)) != groupIDtoArmID.end())
                            UCB1.markForRemoval(groupIDtoArmID[std::make_pair(index1, index2)]);
                        else if (groupIDtoArmID.find(std::make_pair(index2, index1)) != groupIDtoArmID.end())
                            UCB1.markForRemoval(groupIDtoArmID[std::make_pair(index2, index1)]);
                        else
                            std::cout << "Problem when removing right, index " << index1 << " " << index2 << std::endl;
                    }
                }
            }
                for (unsigned index1(0); index1 < m; index1++) {
                    if (groupMembership[index1] == groupMembership[right])
                        groupMembership[index1] = groupMembership[left];
                }
        }
//        for (const auto &index: verticesAdded) {
//            if ( (index==right) or (index==left))
//                continue;
//            if( (left!=index) and (verticesAdded.find(left) == verticesAdded.end()) ){
//                std::cout << left << "," << index << "\t" ;
//                if (groupIDtoArmID.find(std::make_pair(left, index)) != groupIDtoArmID.end())
//                    UCB1.markForRemoval(groupIDtoArmID[std::make_pair(left, index)]);
//                else if (groupIDtoArmID.find(std::make_pair(index, left)) != groupIDtoArmID.end())
//                    UCB1.markForRemoval(groupIDtoArmID[std::make_pair(index, left)]);
//                else
//                    std::cout << "Problem when removing left, index " << left << " " << index << std::endl;
//            }
//
//            if( (right!=index) and (verticesAdded.find(right) == verticesAdded.end()) ){
//                std::cout << right << "," << index << "\t" ;
//                if (groupIDtoArmID.find(std::make_pair(right, index)) != groupIDtoArmID.end())
//                    UCB1.markForRemoval(groupIDtoArmID[std::make_pair(right, index)]);
//                else if (groupIDtoArmID.find(std::make_pair(index, right)) != groupIDtoArmID.end())
//                    UCB1.markForRemoval(groupIDtoArmID[std::make_pair(index, right)]);
//                else
//                    std::cout << "Problem when removing right, index " << right << " " << index << std::endl;
//            }
//        }
        std::cout << std::endl;
        verticesAdded.insert(left);
        verticesAdded.insert(right);
    }

}
