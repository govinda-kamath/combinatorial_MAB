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
#include "../deprecated/UCB.h"
#include "../utilities/INIReader.h"
#include "tenXReader.h"
#include "Knn.h"
#define DEBUG


int main(int argc, char *argv[]){
    std::string nameConfig = argv[1];
    long startIndex(atol(argv[2])); // Start index
    long endIndex(atol(argv[3])); // End index


    INIReader reader(nameConfig);
    if (reader.ParseError() < 0) {
        std::cout << "Can't load "<< nameConfig << std::endl;
        return 1;
    }
    std::string directoryPath = reader.Get("path", "directory", "");
    std::string saveFileFolder =reader.Get("path", "saveFileFolder_10x", "../experiments/knn/tmp10x");
    std::string fileName = reader.Get("path", "h5path", "test_dataset/1M_neurons_matrix_subsampled_2k.h5");
    unsigned numberOfInitialPulls = (unsigned) reader.GetInteger("UCB", "numberOfInitialPulls_knn", 100);
    unsigned sampleSize = (unsigned) reader.GetInteger("UCB", "sampleSize", 32);
    unsigned k = (unsigned) reader.GetInteger("UCB", "k", 5);
    float delta = (float) reader.GetReal("UCB", "delta", 0.1);
    long n = (unsigned) reader.GetInteger("UCB", "n", -1);


    // Loading 10x data shape
    std::vector<unsigned> shapeData =  tenXReader::get10xMatrixSize(fileName);

    // Loading Sparse matrix
    std::cout << "Reading normalised data sparsely" << std::endl;
    std::vector<std::unordered_map<unsigned long, float> > sparseNormalisedDataMatrix(shapeData[1] );
    tenXReader::get10xNormalisedSparseMatrix(fileName, sparseNormalisedDataMatrix);

    //Arms
    std::vector<SparseL1Point> pointsVec;
    utils::unorderedMapToPoints(pointsVec, sparseNormalisedDataMatrix, shapeData[0], n);

    std::cout << "Running "<<k<< "-nn for " << endIndex-startIndex << " points" << std::endl;
    std::cout << "n = " << pointsVec.size() << std::endl;
    std::cout << numberOfInitialPulls << std::endl;

    // Arms and UCB
    std::chrono::system_clock::time_point loopTimeStart = std::chrono::system_clock::now();
    Knn<SparseL1Point> knn(pointsVec, k, numberOfInitialPulls, delta, sampleSize, saveFileFolder);
    std::vector<unsigned long> indices(endIndex-startIndex);
    std::iota(indices.begin(), indices.end(), startIndex);
    std::cout << "Running" <<std::endl;
    knn.run();

    std::chrono::system_clock::time_point loopTimeEnd = std::chrono::system_clock::now();
    std::cout << "Average time (ms) "
              << std::chrono::duration_cast<std::chrono::milliseconds>(loopTimeEnd - loopTimeStart).count()/
                      (endIndex-startIndex) << std::endl;
    return 0;
}

