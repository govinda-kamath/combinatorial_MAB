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
#include "../utilities/INIReader.h"
#include "Knn.h"
#include "HeirarchicalClustering.h"


int main(int argc, char *argv[]) {
//    std::string nameConfig = argv[1];
    std::string nameConfig = "/Users/vivekkumarbagaria/Code/combinatorial_MAB/nominal.ini";

    // Parameters
    INIReader reader(nameConfig);
    if (reader.ParseError() < 0) {
        std::cout << "Can't load " << nameConfig << std::endl;
        return 1;
    }

    // Loading Hyper parameters and data sizes
    std::string directoryPath = reader.Get("path", "h_directory", "");
    std::string graphFilePath = reader.Get("path", "graphFilePathHeirarchical", "h_graph_output");
    std::string saveFilePath = reader.Get("path", "saveFilePathHeirarchical", "h_output_data");
    std::string fileSuffix = reader.Get("path", "suffix", "");
    unsigned numberOfInitialPulls = (unsigned) reader.GetInteger("UCB", "numberOfInitialPulls_knn", 100);
    unsigned sampleSize = (unsigned) reader.GetInteger("UCB", "sampleSize", 32);
    float delta = (float) reader.GetReal("UCB", "delta", 0.1);
    unsigned long n = (unsigned) reader.GetInteger("UCB", "n_h", 100);
    delta = delta / (n);

    // Data to Points (vectors)
    std::vector<std::string> pathsToImages;
    std::vector<SquaredEuclideanPoint> pointsVec;
    utils::getPathToFile(pathsToImages, directoryPath, fileSuffix); // Loads the filepaths into an array
    utils::vectorsToPoints(pointsVec, pathsToImages, n); //Loads the images from the locations in above array to pointsVec

    // Arms and UCB
    std::chrono::system_clock::time_point loopTimeStart = std::chrono::system_clock::now();
    Heirarchical<SquaredEuclideanPoint> heirar(pointsVec, numberOfInitialPulls, delta, sampleSize,
                                               saveFilePath, graphFilePath, n);

    std::cout << "Start!" << std::endl;
    heirar.run();
    std::chrono::system_clock::time_point loopTimeEnd = std::chrono::system_clock::now();

    std::cout << "Average time (ms) "
              << std::chrono::duration_cast<std::chrono::milliseconds>(loopTimeEnd - loopTimeStart).count()/
                 (n) << std::endl;
    heirar.saveAnswers();
}
