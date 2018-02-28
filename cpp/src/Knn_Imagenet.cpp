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
#include "UCB.h"
#include "../utilities/INIReader.h"
#include "Knn.h"
//#define DEBUG


int main(int argc, char *argv[]){
//    std::string nameConfig = argv[1];
//    long startIndex(atol(argv[2])); // Start index
//    long endIndex(atol(argv[3])); // End index
//
//    int tmp(0);

//     For debugging mode in CLion
    std::string nameConfig = "/Users/govinda/Code/combinatorial_MAB/nominal1.ini";
    long startIndex(0); // Start index
    long endIndex(10); // End index

    // Parameters
    INIReader reader(nameConfig);
    if (reader.ParseError() < 0) {
        std::cout << "Can't load "<< nameConfig << std::endl;
        return 1;
    }
    std::string directoryPath = reader.Get("path", "directory", "");
    std::string saveFilePath =reader.Get("path", "saveFilePath", "test.output");
    std::string fileSuffix = reader.Get("path", "suffix", "");
    unsigned numberOfInitialPulls = (unsigned) reader.GetInteger("UCB", "numberOfInitialPulls_knn", 100);
    unsigned k = (unsigned) reader.GetInteger("UCB", "k", 5);
    float delta = (float) reader.GetReal("UCB", "delta", 0.1);

    std::cout << "Running "<<k<< "-nn for " << endIndex-startIndex << " points" << std::endl;
    std::cout << numberOfInitialPulls << std::endl;

    // Data to Points
    std::vector<std::string>  pathsToImages;
    std::vector<SquaredEuclideanPoint> pointsVec;
    utils::getPathToFile(pathsToImages, directoryPath, fileSuffix);
    utils::vectorsToPoints(pointsVec, pathsToImages);

    // Arms and UCB
    std::chrono::system_clock::time_point loopTimeStart = std::chrono::system_clock::now();
    Knn<SquaredEuclideanPoint> knn(pointsVec, k, numberOfInitialPulls, delta );
    std::vector<unsigned long> indices(endIndex-startIndex);
    std::iota(indices.begin(), indices.end(), startIndex);
    std::cout << "Running" << std::endl;
    knn.run(indices);
    std::chrono::system_clock::time_point loopTimeEnd = std::chrono::system_clock::now();
    std::cout << "Saving onw" <<std::endl;
    knn.saveAnswers("./test1.out");

    std::cout << "Average time (ms) "
              << std::chrono::duration_cast<std::chrono::milliseconds>(loopTimeEnd - loopTimeStart).count()/
                      (endIndex-startIndex) << std::endl;

    return 0;
}
