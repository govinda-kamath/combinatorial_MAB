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
#include "Medoid.h"
#include "../utilities/INIReader.h"
#include "utils.h"
#define DEBUG


int main(int argc, char *argv[]){
    std::string nameConfig = argv[1];

//    std::string nameConfig = "nominal.ini"; //For debugging mode in CLion

    INIReader reader(nameConfig);
    if (reader.ParseError() < 0) {
        std::cout << "Can't load "<< nameConfig << std::endl;
        return 1;
    }
    std::string directoryPath = reader.Get("path", "directory", "");
    std::string saveFilePath =reader.Get("path", "saveFilePath", "test.output");
    std::string fileSuffix = reader.Get("path", "suffix", "");
    unsigned numberOfInitialPulls = (unsigned) reader.GetInteger("UCB", "numberOfInitialPulls_medoid", 100);
    float delta = (float) reader.GetReal("UCB", "delta", 0.01);

    std::cout << "Running Medoid" << std::endl;

    // Load data
    std::vector<std::string>  pathsToImages;
    utils::getPathToFile(pathsToImages, directoryPath, fileSuffix);

    // Points
    std::vector<L1Point> pointsVec;
    utils::vectorsToPoints(pointsVec, pathsToImages);

    // Finding the medoids
    unsigned  k = 10;
    Medoid<L1Point> medoid10x(pointsVec, k, numberOfInitialPulls, delta);
    medoid10x.run();
    medoid10x.printOrder();

}
