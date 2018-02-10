//
// Created by Govinda Kamath on 2/8/18.
//

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <thread>
#include <numeric>
#include <mutex>
#include <map>
#include<chrono>
#include <dlib/image_io.h>
#include <dlib/image_transforms.h>
#include "tenXReader.h"

int main(int argc, char *argv[]){
    std::string fileName("/Users/govinda/Code/combinatorial_MAB/test_dataset/1M_neurons_neuron20k.h5");

    std::vector<int> shapeData(2);
    tenXReader::get10xMatrixSize(fileName, shapeData);

    std::chrono::system_clock::time_point dataLoadTimeStart = std::chrono::system_clock::now();
    std::vector<std::vector<int> > denseDataMatrix(shapeData[0], std::vector<int>(shapeData[1]));
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
        for (int j = 0; j < denseDataMatrix.size(); ++j){
            if (denseDataMatrix[j][i] != 0){
                numberNonzeros++;
            }
        }
        std::cout << i << "  " << numberNonzeros << std::endl;
    }


}
