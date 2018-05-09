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


int main(int argc, char *argv[]){

//    std::string nameConfig = argv[1];
//    std::string nameConfig = "/Users/vivekkumarbagaria/Code/combinatorial_MAB/nominal.ini";
//
////     Parameters
//    INIReader reader(nameConfig);
//    if (reader.ParseError() < 0) {
//        std::cout << "Can't load " << nameConfig << std::endl;
//        return 1;
//    }

    // Loading Hyper parameters and data sizes
//    unsigned sampleSize = (unsigned) reader.GetInteger("UCB", "sampleSize", 32);
//    float delta = (float) reader.GetReal("UCB", "delta", 0.1);
//    unsigned long n = (unsigned) reader.GetInteger("UCB", "n_h", 100);
//    delta = delta / (n);

    float delta = 0.1;
    unsigned sampleSize = 100;

    std::random_device rd;
    std::mt19937 g(9);
    std::cout << "Started!" << std::endl;
    std::string str;
    float tmp;

    std::unordered_map<std::pair<unsigned long, unsigned long>, unsigned long, utils::pair_hash> groupIDtoArmID;
    std::vector<std::pair<unsigned long, unsigned long>> armIDtoGroupID;
    //Maintains the valid group points currently and used while creating new arms
    std::unordered_set<unsigned long> verticesAdded;

    long n = 10000;
    long m = 10;
    unsigned long armID = 0;
    std::vector<SquaredEuclideanPoint> dataMatrix;
    std::vector<Arm2DMutualInformation<SquaredEuclideanPoint> > armsVec;

    std::ifstream file("/Users/vivekkumarbagaria/Code/test_dataset/gene_data_MMI.txt", std::ios::binary);
    int ii = 0;
    while (getline(file, str)) {
        std::stringstream ss(str);
        std::cout << ii++ << " ";
        std::vector<float> tmpVec;
        for (long j = 0; j < m; j++){//5000
            ss >> tmp;
            tmpVec.push_back(tmp);
        }
        dataMatrix.push_back(SquaredEuclideanPoint(tmpVec));
        if(ii>1000)
            break;
    }
    std::cout << std::endl;

    armIDtoGroupID.reserve((int)0.5*m*(m-1));
    for(unsigned i(1); i< m ; i ++) {
        for (unsigned j(0); j < i; j++) {
            std::vector<unsigned long> indices = {i, j};
            groupIDtoArmID[std::make_pair(i, j)] = armID;
            armIDtoGroupID.push_back(std::make_pair(i, j));
            std::vector<unsigned long> shuffledRows_(dataMatrix.size());
            std::iota(shuffledRows_.begin(), shuffledRows_.end(), 0);
            std::shuffle(shuffledRows_.begin(), shuffledRows_.end(), g);
            Arm2DMutualInformation<SquaredEuclideanPoint> arm(armID, dataMatrix, indices,shuffledRows_);
            armsVec.push_back(arm);
            std::cout << armID << " " << armIDtoGroupID.size()
                      <<  " " << armIDtoGroupID[armID].first
                       << " " << armIDtoGroupID[armID].second  << std::endl;
            armID++;
        }
    }
    UCBDynamic<Arm2DMutualInformation<SquaredEuclideanPoint> > UCB1(armsVec, delta, 1, 0, sampleSize);
    UCB1.initialise(100);
    for (unsigned long i(0); i < m - 2; i++) {

        //Find the best group points to join
//            ArmHeirarchical<templatePoint> *bestArm;
#ifndef Brute
        UCB1.runUCB(n * 10);
        Arm2DMutualInformation<SquaredEuclideanPoint> bestArm = UCB1.topKArms.back();
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


        std::cout << "Step " << i
                  //                      << "\tTrue. = " << bestArm.trueMean()
                  //                << "\tSecond True. = " << UCB1.topValidArm().trueMean()
                  << "\tEst. = " << bestArm.estimateOfMean
                  << "\tId. 1= " << left
                  << "\tId. 2= " << right
                  << "\tAverage number of Pulls = " << 2*UCB1.globalNumberOfPulls / ((n - 1) * (n - 1))
                  << std::endl;

        // Removing left and right groups from the list of groups (nodes)


        //Removing arms containing either of the above left and right group of points.
        std::cout << "Removing ";
        for (const auto &index: verticesAdded) {
            std::cout << index << " " ;
            if (groupIDtoArmID.find(std::make_pair(left, index)) != groupIDtoArmID.end())
                UCB1.markForRemoval(groupIDtoArmID[std::make_pair(left, index)]);
            else if (groupIDtoArmID.find(std::make_pair(index, left)) != groupIDtoArmID.end())
                UCB1.markForRemoval(groupIDtoArmID[std::make_pair(index, left)]);
            else
                std::cout << "Problem when removing left, index " << left << " " << index << std::endl;

            if (groupIDtoArmID.find(std::make_pair(right, index)) != groupIDtoArmID.end())
                UCB1.markForRemoval(groupIDtoArmID[std::make_pair(right, index)]);
            else if (groupIDtoArmID.find(std::make_pair(index, right)) != groupIDtoArmID.end())
                UCB1.markForRemoval(groupIDtoArmID[std::make_pair(index, right)]);
            else
                std::cout << "Problem when removing right, index " << right << " " << index << std::endl;
        }
        std::cout << std::endl;
        verticesAdded.insert(left);
        verticesAdded.insert(right);
    }

}
