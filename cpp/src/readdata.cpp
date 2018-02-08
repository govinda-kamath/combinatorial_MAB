//
// Created by Govinda Kamath on 2/7/18.
//

//

#include<iostream>
#include<vector>
#include "tenXReader.h"

int main() {

    std::string fileName("/Users/govinda/Code/combinatorial_MAB/test_dataset/1M_neurons_neuron20k.h5");

    std::vector<int> shapeData(2);
    tenXReader::get10xMatrixSize(fileName, shapeData);

    std::vector<std::vector<int> > denseDataMatrix(shapeData[0], std::vector<int>(shapeData[1]));
    std::cout << denseDataMatrix.size() << std::endl;
    std::cout << denseDataMatrix[0].size() << std::endl;
    std::cout << denseDataMatrix[0][100] << std::endl;

    tenXReader::get10xMatrix(fileName, denseDataMatrix);


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
