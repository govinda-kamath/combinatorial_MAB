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

    std::string nameConfig;
    nameConfig = argv[1];
    INIReader reader(nameConfig);
    if (reader.ParseError() < 0) {
        std::cout << "Can't load " << nameConfig << std::endl;
        nameConfig = "/Users/vivekkumarbagaria/Code/combinatorial_MAB/nominal.ini";
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

    for(unsigned i(1); i< m ; i ++) {
        std::chrono::system_clock::time_point sTime = std::chrono::system_clock::now();
        std::vector<unsigned long> indices = {0,i};
        std::vector<unsigned long> shuffledRows_(allPointsVec.size());
        std::iota(shuffledRows_.begin(), shuffledRows_.end(), 0);
        std::shuffle(shuffledRows_.begin(), shuffledRows_.end(), g);
        Arm2DMutualInformation<SquaredEuclideanPoint> arm(armID, allPointsVec, indices, shuffledRows_);
        std::chrono::system_clock::time_point eTime = std::chrono::system_clock::now();
        long long int totalTime = std::chrono::duration_cast<std::chrono::milliseconds>
                (eTime - sTime).count();
        armsVec.push_back(arm);
        std::cout << totalTime << " " << i << std::endl;
    }

    std::cout << "UCB go! for steps " << m-1 << "with arms " << armID << std::endl;
    UCBDynamic<Arm2DMutualInformation<SquaredEuclideanPoint> > UCB1(armsVec, delta, 1, 0, sampleSize);
    UCB1.initialise(numberOfInitialPulls);

    UCB1.runUCB(n * 100);

    Arm2DMutualInformation<SquaredEuclideanPoint> bestArm = UCB1.topKArms.back();
    Arm2DMutualInformation<SquaredEuclideanPoint> secondBestArm = UCB1.arms.top();

}
