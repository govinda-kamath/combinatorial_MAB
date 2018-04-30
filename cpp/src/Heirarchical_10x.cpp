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
#include "tenXReader.h"
#include "../utilities/INIReader.h"
#include "Knn.h"
#include "HeirarchicalClustering.h"


int main(int argc, char *argv[]) {
    std::string nameConfig = argv[1];
    char algo = argv[2][0]; //m - for mab and b for brute


    // Parameters
    INIReader reader(nameConfig);
    if (reader.ParseError() < 0) {
        std::cout << "Can't load " << nameConfig << std::endl;
        return 1;
    }

    // Loading Hyper parameters and data sizes
    std::string directoryPath = reader.Get("path", "h_directory", "");
    std::string fileName = reader.Get("path", "h5path", "test_dataset/1M_neurons_matrix_subsampled_2k.h5");
    std::string graphFilePath = reader.Get("path", "graphFilePathHeirarchical10x", "h_graph_output10x");
    std::string saveFilePath = reader.Get("path", "saveFilePathHeirarchical10x", "h_output_data10x");
    std::string fileSuffix = reader.Get("path", "suffix", "");
    unsigned numberOfInitialPulls = (unsigned) reader.GetInteger("UCB", "numberOfInitialPulls_knn", 100);
    unsigned sampleSize = (unsigned) reader.GetInteger("UCB", "sampleSize", 32);
    float delta = (float) reader.GetReal("UCB", "delta", 0.1);
    unsigned long n = (unsigned) reader.GetInteger("UCB", "n_h", 100);
    delta = delta / (n);

    // Loading 10x data shape
    std::vector<unsigned> shapeData =  tenXReader::get10xMatrixSize(fileName);
    // Loading Sparse matrix
    std::cout << "Reading normalised data sparsely" << std::endl;
    std::vector<std::unordered_map<unsigned long, float> > sparseNormalisedDataMatrix(shapeData[1] );
    tenXReader::get10xNormalisedSparseMatrix(fileName, sparseNormalisedDataMatrix);
    //Arms
    std::vector<SparseL1Point> pointsVec;
    utils::unorderedMapToPoints(pointsVec, sparseNormalisedDataMatrix, shapeData[0]);

    //UCB
    std::chrono::system_clock::time_point loopTimeStart = std::chrono::system_clock::now();
    Heirarchical<SparseL1Point> heirar(pointsVec, numberOfInitialPulls, delta, sampleSize,
                                               saveFilePath, graphFilePath,  algo, n);

    std::cout << "Start!" << std::endl;
    heirar.run();
    std::chrono::system_clock::time_point loopTimeEnd = std::chrono::system_clock::now();

    std::cout << "Average time (ms) "
              << std::chrono::duration_cast<std::chrono::milliseconds>(loopTimeEnd - loopTimeStart).count()/
                 (n) << std::endl;
    heirar.saveAnswers();
}
