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
    std::string nameConfig = "/Users/vivekkumarbagaria/Code/combinatorial_MAB/nominal.ini";
//
////     Parameters
    INIReader reader(nameConfig);
    if (reader.ParseError() < 0) {
        std::cout << "Can't load " << nameConfig << std::endl;
        return 1;
    }

    std::random_device rd;
    std::mt19937 g(9);
    std::normal_distribution<> dist(0, 0.0005);
    // Loading Hyper parameters and data sizes
    unsigned sampleSize = (unsigned) reader.GetInteger("UCB", "sampleSize", 32);
    float delta = (float) reader.GetReal("UCB", "delta", 0.1);
    unsigned numberOfInitialPulls = (unsigned) reader.GetInteger("UCB", "numberOfInitialPulls_knn", 100);
    unsigned long m = (unsigned) reader.GetInteger("UCB", "m", 100);
    unsigned long n = (unsigned) reader.GetInteger("UCB", "n", 100);
//    delta = delta / (n);

//    float delta = 0.1;
//    unsigned sampleSize = 100;


    std::cout << "Started!" << std::endl;
    std::string str;
    float tmp;

    std::unordered_map<std::pair<unsigned long, unsigned long>, unsigned long, utils::pair_hash> groupIDtoArmID;
    std::vector<std::pair<unsigned long, unsigned long>> armIDtoGroupID;
    //Maintains the valid group points currently and used while creating new arms
    std::unordered_set<unsigned long> verticesAdded;

//    long n = 10000;
    unsigned long armID = 0;
    std::vector<SquaredEuclideanPoint> dataMatrix;
    std::vector<Arm1DEntropy<SquaredEuclideanPoint> > armsVec;

    std::ifstream file("/Users/vivekkumarbagaria/Code/test_dataset/gene_data_MMI.txt", std::ios::binary);
    int ii = 0;
    while (getline(file, str)) {
        std::stringstream ss(str);
        if(ii%100==0)
            std::cout << ii << std::endl;
        std::vector<float> tmpVec;
        tmpVec.reserve(m);
        for (long j = 0; j < m; j++){//5000
            ss >> tmp;
            if(tmp!=0)
                tmp += dist(g);
            tmpVec.push_back(tmp);
        }
        dataMatrix.push_back(SquaredEuclideanPoint(tmpVec));
        if(ii>n)
            break;
        ii++;
    }

    armIDtoGroupID.reserve((int)0.5*m*(m-1));
    for(unsigned i(0); i< m ; i ++) {
        std::vector<unsigned long> shuffledRows_(dataMatrix.size());
        std::iota(shuffledRows_.begin(), shuffledRows_.end(), 0);
        std::shuffle(shuffledRows_.begin(), shuffledRows_.end(), g);
        Arm1DEntropy<SquaredEuclideanPoint> arm(armID, dataMatrix, i, shuffledRows_);
        armsVec.push_back(arm);
        armID++;
    }

    std::cout << "UCB go! for steps " << m-1 << std::endl;
    UCBDynamic<Arm1DEntropy<SquaredEuclideanPoint> > UCB1(armsVec, delta, 1, 0, sampleSize);
    UCB1.initialise(numberOfInitialPulls);

    for (unsigned long i(0); i < 0; i++) {
        UCB1.runUCB(n * 10);
        Arm1DEntropy<SquaredEuclideanPoint> bestArm = UCB1.topKArms.back();
        Arm1DEntropy<SquaredEuclideanPoint> secondBestArm = UCB1.arms.top();
    }

}
