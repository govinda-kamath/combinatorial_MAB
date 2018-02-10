//
// Created by Govinda Kamath on 2/7/18.
//

//

#include<iostream>
#include<vector>
#include<chrono>
#include "tenXReader.h"

int main() {

    std::string fileName("/Users/govinda/Code/combinatorial_MAB/test_dataset/1M_neurons_neuron20k.h5");

    std::vector<int> shapeData(2);
    tenXReader::get10xMatrixSize(fileName, shapeData);

    std::chrono::system_clock::time_point dataLoadTimeStart = std::chrono::system_clock::now();
    std::vector<std::vector<int> > denseDataMatrix(shapeData[1], std::vector<int>(shapeData[0]));
    std::cout << denseDataMatrix.size() << std::endl;
    std::cout << denseDataMatrix[0].size() << std::endl;
    std::cout << denseDataMatrix[0][100] << std::endl;

    tenXReader::get10xMatrix(fileName, denseDataMatrix);
    std::chrono::system_clock::time_point dataLoadTimeEnd = std::chrono::system_clock::now();
    std::cout << "Time to read dense vector (ms) "
              << std::chrono::duration_cast<std::chrono::milliseconds>(dataLoadTimeEnd - dataLoadTimeStart).count()
              << std::endl;

    for (int i = 0; i < 10; ++i) {
        int numberNonzeros(0);
        for (int j = 0; j < denseDataMatrix[0].size(); ++j){
            if (denseDataMatrix[i][j] != 0){
                numberNonzeros++;
            }
        }
        std::cout << i << "  " << numberNonzeros << std::endl;
    }

    std::cout << "Reading data normalisec " << std::endl;
    dataLoadTimeStart = std::chrono::system_clock::now();
    std::vector<std::vector<float> > denseNormalisedDataMatrix(shapeData[1], std::vector<float>(shapeData[0]));
    std::cout << denseNormalisedDataMatrix.size() << std::endl;
    std::cout << denseNormalisedDataMatrix[0].size() << std::endl;
    std::cout << denseNormalisedDataMatrix[0][100] << std::endl;

    tenXReader::get10xNormalisedDenseMatrix(fileName, denseNormalisedDataMatrix);
    dataLoadTimeEnd = std::chrono::system_clock::now();
    std::cout << "Time to read normalised dense vector (ms) "
              << std::chrono::duration_cast<std::chrono::milliseconds>(dataLoadTimeEnd - dataLoadTimeStart).count()
              << std::endl;

    for (int i = 0; i < 10; ++i) {
        int numberNonzeros(0);
        for (int j = 0; j < denseNormalisedDataMatrix[0].size(); ++j){
            if (denseNormalisedDataMatrix[i][j] != 0){
                numberNonzeros++;
            }
        }
        std::cout << i << "  " << numberNonzeros << std::endl;
    }


    std::cout << "Reading data sparsely " << std::endl;
    dataLoadTimeStart = std::chrono::system_clock::now();
    std::vector<std::unordered_map<unsigned long, int> > sparseDataMatrix(shapeData[1] );
    tenXReader::get10xSparseMatrix(fileName, sparseDataMatrix);
    dataLoadTimeEnd = std::chrono::system_clock::now();
    std::cout << "Time to read sparse vector (ms) "
              << std::chrono::duration_cast<std::chrono::milliseconds>(dataLoadTimeEnd - dataLoadTimeStart).count()
              << std::endl;
    for (int i = 0; i < 10; ++i) {
        std::cout << i << "  " << sparseDataMatrix[i].size() << std::endl;
    }

    std::cout << "Reading normalised data sparsely " << std::endl;
    dataLoadTimeStart = std::chrono::system_clock::now();
    std::vector<std::unordered_map<unsigned long, float> > sparseNormalisedDataMatrix(shapeData[1] );
    tenXReader::get10xNormalisedMatrix(fileName, sparseNormalisedDataMatrix);
    dataLoadTimeEnd = std::chrono::system_clock::now();
    std::cout << "Time to read normalised sparse vector (ms) "
              << std::chrono::duration_cast<std::chrono::milliseconds>(dataLoadTimeEnd - dataLoadTimeStart).count()
              << std::endl;
    for (int i = 0; i < 10; ++i) {
        std::cout << i << "  " << sparseDataMatrix[i].size() << std::endl;
    }

}
