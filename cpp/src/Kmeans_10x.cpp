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
#include "INIReader.h"
#include "kmeans.h"
#include "tenXReader.h"
#define DEBUG

int main(int argc, char *argv[]){
    std::string nameConfig = argv[1];
//    long startIndex(atol(argv[2])); // Start index
//    long endIndex(atol(argv[3])); // End index


//     For debugging mode in CLion
//    std::string nameConfig = "nominal.ini";

    INIReader reader(nameConfig);
    if (reader.ParseError() < 0) {
        std::cout << "Can't load "<< nameConfig << std::endl;
        return 1;
    }
    std::string directoryPath = reader.Get("path", "directory", "");
    std::string saveFileFolder =reader.Get("path", "saveFileFolderkmean10x", "../experiments/kmeans/10x/");
    std::string fileName = reader.Get("path", "h5path", "test_dataset/1M_neurons_matrix_subsampled_2k.h5");
    unsigned numberOfInitialPulls = (unsigned) reader.GetInteger("UCB", "numberOfInitialPulls_knn", 100);
    unsigned sampleSize = (unsigned) reader.GetInteger("UCB", "sampleSize", 32);
    unsigned k = (unsigned) reader.GetInteger("UCB", "k", 5);
    float delta = (float) reader.GetReal("UCB", "delta", 0.1);
    long n = (unsigned) reader.GetInteger("UCB", "n", -1);
    delta = delta / (k);


    // Loading 10x data shape
    std::vector<unsigned> shapeData =  tenXReader::get10xMatrixSize(fileName);

    // Loading Sparse matrix
    std::cout << "Reading normalised data sparsely" << std::endl;
    std::vector<std::unordered_map<unsigned long, float> > sparseNormalisedDataMatrix(shapeData[1] );
    tenXReader::get10xNormalisedSparseMatrix(fileName, sparseNormalisedDataMatrix);

    //Arms
    std::vector<SparseL1Point> allPointsVec;
    utils::unorderedMapToPoints(allPointsVec, sparseNormalisedDataMatrix, shapeData[0], n);
    
    std::vector<SparseL1Point> centersVec;
    std::vector<SparseL1Point> pointsVec;
    if (n>allPointsVec.size())
        n = allPointsVec.size();
    for(unsigned long i(0); i < n ; i++)
        pointsVec.push_back(allPointsVec[i]);
    //Choosing random points as my cluster centers
    for(unsigned long i(0); i < k ; i++){
        unsigned long index = std::rand()%n;
        centersVec.push_back(allPointsVec[index]);
    }

    //Running one step of kmeans, the maximization step.
    std::chrono::system_clock::time_point loopTimeStart = std::chrono::system_clock::now();
    Kmeans<SparseL1Point> kMeans( pointsVec, centersVec, numberOfInitialPulls, delta, sampleSize, saveFileFolder);
    kMeans.maximization();
    std::chrono::system_clock::time_point loopTimeEnd = std::chrono::system_clock::now();

    std::cout << "Average time (ms) "
              << std::chrono::duration_cast<std::chrono::milliseconds>(loopTimeEnd - loopTimeStart).count()/
                 (k) << std::endl;
    kMeans.saveAnswers();

}


