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
#define DEBUG

int main(int argc, char *argv[]){
//    std::string nameConfig = argv[1];
//    long startIndex(atol(argv[2])); // Start index
//    long endIndex(atol(argv[3])); // End index


//     For debugging mode in CLion
    std::string nameConfig = "nominal.ini";
    long startIndex(0); // Start index
    long endIndex(100); // End index

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


    std::cout << "Running "<<k<< "-means for " << endIndex-startIndex << " points" << std::endl;
    std::cout << numberOfInitialPulls << std::endl;

    // Data
    std::vector<std::string>  pathsToImages;
    utils::getPathToFile(pathsToImages, directoryPath, fileSuffix);

    //Points
    std::vector<SquaredEuclideanPoint> initialPointsVec;
    utils::vectorsToPoints(initialPointsVec, pathsToImages);

    // Arms and UCB
    std::vector<SquaredEuclideanPoint> centersVec;
    std::vector<SquaredEuclideanPoint> pointsVec;

    for(unsigned long i(0); i < k ; i++){
        unsigned long index = std::rand()%initialPointsVec.size();
        centersVec.push_back(initialPointsVec[index]);
    }
    for(unsigned long i(k); i < initialPointsVec.size() ; i++)
        pointsVec.push_back(initialPointsVec[i]);

    std::chrono::system_clock::time_point loopTimeStart = std::chrono::system_clock::now();
    Kmeans<SquaredEuclideanPoint> kMeans( pointsVec, centersVec, numberOfInitialPulls, delta);
    std::vector<unsigned long> indices(endIndex-startIndex);
    std::iota(indices.begin(), indices.end(), startIndex);
//    kmeans.run(indices);
    kMeans.run();


    std::chrono::system_clock::time_point loopTimeEnd = std::chrono::system_clock::now();
    std::cout << "Average time (ms) "
              << std::chrono::duration_cast<std::chrono::milliseconds>(loopTimeEnd - loopTimeStart).count()/
                 (endIndex-startIndex) << std::endl;

    loopTimeStart = std::chrono::system_clock::now();
    for (unsigned  long i(0); i << centersVec.size(); i++) {
        for (unsigned long j(0); j << pointsVec.size(); j++) {
            float die = centersVec[i].distance(pointsVec[j]);
        }
    }
    loopTimeEnd = std::chrono::system_clock::now();

    std::cout << "Average time Naive (ms) "
              << std::chrono::duration_cast<std::chrono::milliseconds>(loopTimeEnd - loopTimeStart).count()/
                 (endIndex-startIndex) << std::endl;
     return 0;
}