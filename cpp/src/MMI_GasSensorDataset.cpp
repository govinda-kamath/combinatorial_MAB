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
    std::chrono::system_clock::time_point timeS = std::chrono::system_clock::now();
    std::string nameConfig;
    nameConfig = argv[1];
    long m(atol(argv[2]));
    long n(atol(argv[3]));
    long rrr(atol(argv[4]));
    std::srand(std::time(nullptr));
//    nameConfig = "/Users/vivekkumarbagaria/Code/combinatorial_MAB/nominal.ini";

    INIReader reader(nameConfig);
    if (reader.ParseError() < 0) {
        std::cout << "Can't load " << nameConfig << std::endl;
    }

//    Loading Hyper parameters and data sizes
    unsigned sampleSize = (unsigned) reader.GetInteger("UCB", "sampleSize", 32);
    float delta = (float) reader.GetReal("UCB", "delta", 0.1);
    unsigned numberOfInitialPulls = (unsigned) reader.GetInteger("UCB", "numberOfInitialPulls_knn", 100);
//    unsigned long m = (unsigned) reader.GetInteger("UCB", "m", 100);
//    unsigned long n = (unsigned) reader.GetInteger("UCB", "n", 100);
    std::string fileName = reader.Get("path", "h5path", "../test_dataset/1M_neurons_matrix_subsampled_20k_5kgenes.h5");
    std::string directoryPath = reader.Get("path", "gas_directory", "");
    std::string saveFilePath = reader.Get("path", "saveFileFolderGas", "../experiments/MI/Gas/");


    std::cout << "Started!" << std::endl;
    std::string str;
    float tmp;
    std::vector<Arm2DMutualInformation<SquaredEuclideanPoint> > armsVec;
//
//    // Loading 10x data shape
//    std::vector<unsigned> shapeData =  tenXReader::get10xMatrixSize(fileName);
//
//    // Loading Sparse matrix
//    std::cout << "Reading normalised data" << std::endl;
////    std::vector<std::unordered_map<unsigned long, float> > sparseNormalisedDataMatrix(shapeData[1] );
//    unsigned cells(100000), genes(5000);
//    std::vector<std::vector<float> > DataMatrix(cells);
//    for(unsigned i(0); i< cells ; i ++)
//        DataMatrix[i] = std::vector<float>(genes);
//
//    tenXReader::get10xMatrix(fileName, DataMatrix);
//
//    //Arms
//    std::vector<SquaredEuclideanPoint> allPointsVec;
//
//    std::cout << "Reading normalised data sparsely part 2" << std::endl;
//    utils::vectorsToPoints(allPointsVec, DataMatrix,  n, m);


    std::vector<SquaredEuclideanPoint> allPointsVec;
    std::ifstream file(directoryPath, std::ios::binary);

    int i = 0;
    while (getline(file, str)) {
        std::stringstream ss(str);
        std::vector<float> tmpVec;
        for (long j = 0; j < m; j++){//5000
            ss >> tmp;
//            tmp += dist(g);
            tmpVec.push_back(tmp);
        }
        allPointsVec.push_back(SquaredEuclideanPoint(tmpVec));
        i++;
        if(i>n)
            break;
    }
    std::string sFilePath = saveFilePath+std::to_string(std::rand()%1000)+"n_"+std::to_string(i)+"_d_"+std::to_string(m);

    unsigned mainCol = 0;
    std::mt19937 g(std::time(nullptr));
    std::vector<unsigned long> shuffledRows_(allPointsVec.size());
    std::iota(shuffledRows_.begin(), shuffledRows_.end(), 0);
    std::shuffle(shuffledRows_.begin(), shuffledRows_.end(), g);

    for(unsigned i(0); i< m ; i ++) {
        if(i==mainCol)
            continue;
        std::vector<unsigned long> indices = {mainCol, i};
        Arm2DMutualInformation<SquaredEuclideanPoint> arm(i, allPointsVec, indices, shuffledRows_);
//        std::chrono::system_clock::time_point eTime = std::chrono::system_clock::now();
//        long long int totalTime = std::chrono::duration_cast<std::chrono::milliseconds>
//                (eTime - sTime).count();
        armsVec.push_back(arm);
//        std::cout << totalTime << " " << i << std::endl;
    }


    std::cout << "UCB go! for steps " << m-1 << " with arms " << m << std::endl;
    std::chrono::system_clock::time_point timeE = std::chrono::system_clock::now();
    long long int tTime = std::chrono::duration_cast<std::chrono::nanoseconds>(timeE - timeS).count();
//    std::srand(tTime);
    UCBDynamic<Arm2DMutualInformation<SquaredEuclideanPoint> > UCB1(armsVec, delta, 1, 0, sampleSize);
    std::cout<<"Going for initialization" << std::endl;
    std::chrono::system_clock::time_point timeStart = std::chrono::system_clock::now();

    UCB1.initialise(numberOfInitialPulls);
    std::cout<<"Running UCB" << std::endl;
    UCB1.runUCB(n * 100);
    std::chrono::system_clock::time_point timeEnd = std::chrono::system_clock::now();

    Arm2DMutualInformation<SquaredEuclideanPoint> bestArm = UCB1.topKArms.back();
    Arm2DMutualInformation<SquaredEuclideanPoint> secondBestArm = UCB1.arms.top();

    auto id = bestArm.id;
    auto id2 = secondBestArm.id;


    std::cout<< "Saving in " << sFilePath << std::endl;
    std::ofstream saveFile(sFilePath, std::ofstream::out | std::ofstream::trunc);
    long long int totalTime = std::chrono::duration_cast<std::chrono::milliseconds>
            (timeEnd - timeStart).count();
    saveFile   << "n," << allPointsVec.size() << std::endl;
    saveFile   << "m," << m << std::endl;
    saveFile   << "id," << id << std::endl;
    saveFile   << "id2," << id2 << std::endl;
    saveFile   << "Average number of Pulls," << UCB1.globalNumberOfPulls/UCB1.arms.size()<< std::endl;
    saveFile   << "Gain," << (float) allPointsVec.size()/(UCB1.globalNumberOfPulls/UCB1.arms.size())<< std::endl;
    saveFile   << "secondBestId," << id2 << std::endl;
    saveFile   << "TotalTime," << totalTime << std::endl;
    saveFile   << "GlobalnumberofPulls," << UCB1.globalNumberOfPulls << std::endl;
    saveFile   << "GlobalSigma," << UCB1.globalSigma  << std::endl;
    saveFile   << "Shuffled rows"
            << "," << shuffledRows_[0]
            << "," << shuffledRows_[1]
            << "," << shuffledRows_[2]
            << "," << shuffledRows_[3]
            << std::endl;



    std::cout << "Total Time = " << totalTime << "(ms)"
              << "\nGlobal number of Pulls = " << UCB1.globalNumberOfPulls
              << "\nAverage number of Pulls = " << UCB1.globalNumberOfPulls/UCB1.arms.size()
              << "\nGain, = " << allPointsVec.size()/(UCB1.globalNumberOfPulls/UCB1.arms.size())
              << "\nm," << m
              << "\nn," << allPointsVec.size()
              << "\nGlobal Sigma= " << UCB1.globalSigma
              << std::endl;
    std::cout   << "Shuffled rows"
               << "," << shuffledRows_[0]
               << "," << shuffledRows_[1]
               << "," << shuffledRows_[2]
               << "," << shuffledRows_[3]
            << std::endl;

}
