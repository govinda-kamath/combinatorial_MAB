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
#include "INIReader.h"
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
    unsigned numberOfInitialPulls = (unsigned) reader.GetInteger("UCB", "numberOfInitialPulls", 100);
    float delta = (float) reader.GetReal("UCB", "delta", 0.01);

    std::cout << "Running Medoid" << std::endl;
    std::cout << numberOfInitialPulls << std::endl;

    std::vector<std::string>  pathsToImages;
    utils::getPathToFile(pathsToImages, directoryPath, fileSuffix);

    std::vector<L1Point> pointsVec;
    std::vector<ArmMedoid<L1Point> > armsVec;

    // Points
    for  (unsigned long i(0); i < pathsToImages.size(); i++) {
        std::vector<float> tmpVec;
        utils::readImageAsVector(pathsToImages[i],tmpVec);
        L1Point tmpPoint(tmpVec);
        pointsVec.push_back(tmpPoint);
        if (i%10000 == 9999){
            std::cout << i+1 << " points read." << std::endl;
        }
    }

    //Arms
    for (unsigned i(0); i < pointsVec.size(); i++) {
        ArmMedoid<L1Point> tmpArm(i, pointsVec[i], pointsVec);
        armsVec.push_back(tmpArm);
    }

    std::chrono::system_clock::time_point loopTimeStart = std::chrono::system_clock::now();

    //Running UCB
    unsigned k = 10;
    UCB<ArmMedoid<L1Point> > UCB1(armsVec, delta, k);
    UCB1.initialise(numberOfInitialPulls);
    UCB1.runUCB(200*pointsVec.size());

    //Result
    std::chrono::system_clock::time_point loopTimeEnd = std::chrono::system_clock::now();
    std::cout << "Average time (ms) "
    << std::chrono::duration_cast<std::chrono::milliseconds>(loopTimeEnd - loopTimeStart).count() << std::endl;
    std::cout<<" Top Medoids are";
    for (unsigned i = 0; i < k; i++) {
        std::cout <<  " " << UCB1.topKArms[i].id;
    }
    std::cout<<std::endl;
    return 0;
}
