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
    std::string nameConfig = argv[1];
//    long startIndex(atol(argv[2])); // Start index
//    long endIndex(atol(argv[3])); // End index


//     For debugging mode in CLion
//    std::string nameConfig = "nominal.ini";
    long startIndex(0); // Start index
    long endIndex(100); // End index

    INIReader reader(nameConfig);
    if (reader.ParseError() < 0) {
        std::cout << "Can't load "<< nameConfig << std::endl;
        return 1;
    }
    std::string directoryPath = reader.Get("path", "directory", "");
    std::string saveFileFolder =reader.Get("path", "saveFileFolderkmean", "../experiments/kmeans/imagement/");
    std::string fileSuffix = reader.Get("path", "suffix", "");
    unsigned numberOfInitialPulls = (unsigned) reader.GetInteger("UCB", "numberOfInitialPulls_knn", 100);
    unsigned k = (unsigned) reader.GetInteger("UCB", "k", 50);
    float delta = (float) reader.GetReal("UCB", "delta", 0.1);
    unsigned sampleSize = (unsigned) reader.GetInteger("UCB", "sampleSize", 32);
    long n = (unsigned) reader.GetInteger("UCB", "n", -1);
    delta = delta / (k);


    std::cout << "Running "<<k<< "-means for " << endIndex-startIndex << " points" << std::endl;
    std::cout << numberOfInitialPulls << std::endl;

    // Data
    std::vector<std::string>  pathsToImages;
    utils::getPathToFile(pathsToImages, directoryPath, fileSuffix);

    //Points
    std::vector<SquaredEuclideanPoint> allPointsVec;
    utils::vectorsToPoints(allPointsVec, pathsToImages, n);
    std::vector<SquaredEuclideanPoint> centersVec;
    std::vector<SquaredEuclideanPoint> pointsVec;
    for(unsigned long i(0); i < allPointsVec.size() ; i++)
        pointsVec.push_back(allPointsVec[i]);
    //Choosing random points as my cluster centers
    for(unsigned long i(0); i < k ; i++){
        unsigned long index = std::rand()%allPointsVec.size();
        centersVec.push_back(allPointsVec[index]);
    }

    //Running one step of kmeans, the maximization step.
    std::chrono::system_clock::time_point loopTimeStart = std::chrono::system_clock::now();
    Kmeans<SquaredEuclideanPoint> kMeans( pointsVec, centersVec, numberOfInitialPulls, delta, sampleSize, saveFileFolder);
    kMeans.maximization();
    std::chrono::system_clock::time_point loopTimeEnd = std::chrono::system_clock::now();

    std::cout << "Average time (ms) "
              << std::chrono::duration_cast<std::chrono::milliseconds>(loopTimeEnd - loopTimeStart).count()/
                 (endIndex-startIndex) << std::endl;
    kMeans.saveAnswers();

    // Print the clusters
//    for (unsigned i(0); i < k; i ++){
//        std::cout<< "Cluster " << i  << std::endl;
//        for(unsigned long j(0); j < kMeans.clusters[i].size(); j ++){
//            std::cout<< i << " " << kMeans.clusters[i][j];
//        }
//        std::cout<< std::endl;
//    }

}


