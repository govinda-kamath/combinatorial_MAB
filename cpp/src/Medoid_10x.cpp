//
// Created by Govinda Kamath on 2/8/18.
//

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <thread>
#include <numeric>
#include <mutex>
#include <map>
#include <chrono>
#include <dlib/image_io.h>
#include <dlib/image_transforms.h>
#include "tenXReader.h"
#include "Points.h"
#include "Arms.h"
#include "UCB.h"
#include "INIReader.h"
#include "utils.h"
#include "INIReader.h"
#include "Medoid.h"

int main(int argc, char *argv[]){


    std::string nameConfig = argv[1];
    std::srand(std::time(nullptr));

    INIReader reader(nameConfig);
    if (reader.ParseError() < 0) {
        std::cout << "Can't load "<< nameConfig << std::endl;
        return 1;
    }

    // Reading parameters
    std::string saveFilePath =reader.Get("path", "saveFilePath", "test.output");
    unsigned numberOfInitialPulls = (unsigned) reader.GetInteger("UCB", "numberOfInitialPulls_medoid", 100);
    float delta = (float) reader.GetReal("UCB", "delta", 0.01);
    std::string fileName = reader.Get("path", "h5path", "test_dataset/1M_neurons_matrix_subsampled_2k.h5");
    unsigned k = 3;

    // Loading 10x data shape
    std::vector<unsigned> shapeData =  tenXReader::get10xMatrixSize(fileName);

    // Loading Dense matrix
//    std::vector<std::vector<float> > denseDataMatrix(shapeData[1], std::vector<float>(shapeData[0]));
//    std::cout << "The dense matrix is " << denseDataMatrix.size() << " x " << denseDataMatrix[0].size() << std::endl;
//    tenXReader::get10xNormalisedDenseMatrix(fileName, denseDataMatrix);

//     Loading Sparse matrix
    std::cout << "Reading normalised data sparsely " << std::endl;
    std::vector<std::unordered_map<unsigned long, float> > sparseNormalisedDataMatrix(shapeData[1] );
    tenXReader::get10xNormalisedSparseMatrix(fileName, sparseNormalisedDataMatrix);


    //Arms
    std::vector<SparseL1Point> pointsVec;
    utils::unorderedMapToPoints(pointsVec, sparseNormalisedDataMatrix, shapeData[0]);

//    std::vector<L1Point> pointsVec;
//    utils::vectorsToPoints(pointsVec, denseDataMatrix);
//    std::vector<ArmMedoid<L1Point> > armsVec;

    Medoid<SparseL1Point> medoid10x = Medoid( pointsVec, k, numberOfInitialPulls, delta );
    medoid10x.run();
    medoid10x.printOrder();
}
