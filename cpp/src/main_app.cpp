#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <thread>
#include <map>
#include <dlib/image_io.h>
#include <dlib/image_transforms.h>
#include <glob.h>
#include "Points.h"
#include "Arms.h"
#include "UCB.h"
#include "INIReader.h"
#define DEBUG


//Reads the protected container object of a priority queue to give access to its elements
//A beautiful hack picked up from https://stackoverflow.com/a/12886393/3377314
template <class T, class S, class C>
S& Container(std::priority_queue<T, S, C>& q) {
    struct HackedQueue : private std::priority_queue<T, S, C> {
        static S& Container(std::priority_queue<T, S, C>& q) {
            return q.*&HackedQueue::c;
        }
    };
    return HackedQueue::Container(q);
}



void readImageAsVector (std::string filePath, std::vector<float> &imageVec) {

    dlib::array2d <dlib::rgb_pixel> imageRGB;
    dlib::load_image(imageRGB, filePath.c_str());
    unsigned  numColumns(imageRGB.nc()), numRows(imageRGB.nr());
    unsigned numPixels(numColumns*numRows);
    unsigned vecLength(numPixels*3);

//    std::cout << numColumns <<"\t" << numRows <<
//              "\t" << numPixels << "\t" << vecLength << std::endl;

    if (imageVec.size() != vecLength){

//        std::cout << "initialising" << std::endl;
        imageVec.clear();
        imageVec.reserve(vecLength);

        for (unsigned i(0); i < numRows; i++){
            for (unsigned j(0); j < numColumns; j++) {
                imageVec.push_back((float) imageRGB[i][j].red);
            }
        }



        for (unsigned i(0); i < numRows; i++){
            for (unsigned j(0); j < numColumns; j++) {
                imageVec.push_back((float) imageRGB[i][j].blue);
            }
        }

        for (unsigned i(0); i < numRows; i++){
            for (unsigned j(0); j < numColumns; j++) {
                imageVec.push_back((float) imageRGB[i][j].green);
            }
        }

    } else {

        for (unsigned i(0); i < numRows; i++) {
            for (unsigned j(0); j < numColumns; j++) {
                imageVec[i * numRows + j] = (float) imageRGB[i][j].red;
            }
        }

        for (unsigned i(0); i < numRows; i++) {
            for (unsigned j(0); j < numColumns; j++) {
                imageVec[numPixels + i * numRows + j] = (float) imageRGB[i][j].blue;
            }
        }

        for (unsigned i(0); i < numRows; i++) {
            for (unsigned j(0); j < numColumns; j++) {
                imageVec[2 * numPixels + i * numRows + j] = (float) imageRGB[i][j].green;
            }
        }
    }
}


void singleRun(std::vector<SquaredEuclideanPoint> &pointsVec, unsigned threadNumber, unsigned long mainPointIndexStart,
               unsigned long mainPointIndexEnd, int numberOfInitialPulls, float delta){

    std::cout << "Running thread "<< threadNumber << std::endl;
    for (unsigned long index = mainPointIndexStart; index<mainPointIndexEnd; index++){
        std::chrono::system_clock::time_point loopTimeStart = std::chrono::system_clock::now();
//        std::this_thread::sleep_for(std::chrono::seconds(1*(threadNumber+1)));


        std::vector<ArmKNN<SquaredEuclideanPoint> > armsVec;
        for (unsigned i(0); i < pointsVec.size(); i++) {
            if (i == index)
                continue;
            ArmKNN<SquaredEuclideanPoint> tmpArm(i, pointsVec[i], pointsVec[index]);
            armsVec.push_back(tmpArm);
        }
        std::chrono::system_clock::time_point loopTimeMid = std::chrono::system_clock::now();

        UCB<ArmKNN<SquaredEuclideanPoint> > UCB1(armsVec, delta);
        UCB1.initialise(numberOfInitialPulls);
        UCB1.runUCB(1000000);


        std::chrono::system_clock::time_point loopTimeEnd = std::chrono::system_clock::now();
        if (index%5==0) {
            std::cout << "Ran thread " << threadNumber << ".Index " << index << " Total Time taken = " <<
                      std::chrono::duration_cast<std::chrono::milliseconds>(loopTimeEnd - loopTimeStart).count()
                      << std::endl;
            std::cout << "Ran thread " << threadNumber << ".Index " << index << " Mid Time taken = " <<
                      std::chrono::duration_cast<std::chrono::milliseconds>(loopTimeMid - loopTimeStart).count()
                      << std::endl;
        }
    }

}

int main(int argc, char *argv[]){
    std::cout << "We have entered " << argc
              << " arguments." << std::endl;

    std::vector<SquaredEuclideanPoint> pointsVec;

    std::string nameConfig;
    if (argc==2){
        std::cout<<argc<<std::endl;
        nameConfig = argv[1];
    }
    else{
        nameConfig = "nominal.ini";
    }
    INIReader reader(nameConfig);

    if (reader.ParseError() < 0) {
        std::cout << "Can't load "<< nameConfig << std::endl;
        return 1;
    }

    std::string directoryPath = reader.Get("path", "directory", "");
    std::string fileSuffix = reader.Get("path", "suffix", "");

    int numberOfInitialPulls = (int) reader.GetInteger("UCB", "numberOfInitialPulls", 100);
    float delta = (float) reader.GetReal("UCB", "delta", 0.1);
    int numCores = (int) reader.GetReal("UCB", "numCores", 1);
    unsigned long  maxNumberOfPoints = (unsigned long) reader.GetReal("UCB", "noPoints", 1);

    std::cout << "Running K-nn for " << maxNumberOfPoints << " points using ";
    std::cout << "Number of cores = " << numCores<<std::endl;
    std::cout << numberOfInitialPulls << std::endl;
    std::cout << delta << std::endl;
    std::cout << directoryPath << std::endl;
    std::cout << fileSuffix << std::endl;
    glob_t glob_result;
    std::vector<float> tmpVec;

    unsigned long fileNumber(0);
    std::string searchName;
//
    std::string filePrefix = "/*";
    searchName = directoryPath + filePrefix + fileSuffix;
    std::cout << searchName << std::endl;

    clock_t timeRead = clock();
    glob(searchName.c_str(), GLOB_TILDE, NULL, &glob_result);


    std::cout << "Number of files " << glob_result.gl_pathc << std::endl;
    std::cout << glob_result.gl_pathv[0] << std::endl;

    std::vector<std::string> pathsToImages;
    pathsToImages.reserve(glob_result.gl_pathc);
    for (unsigned long i = 0; i < glob_result.gl_pathc; i ++){
        pathsToImages.push_back(glob_result.gl_pathv[i]);
    }

    std::sort(pathsToImages.begin(), pathsToImages.end());
//    std::map<int,std::string> indexToPath;
//
//
//    for (unsigned long i = 0; i < 100; i ++){
//        indexToPath[i] = pathsToImages[i] ;
//    }

//    while (glob_result.gl_pathc != 0){
//            std::cout << std::string(glob_result.gl_pathv[0]) << std::endl;
//
//        pathsToImages.push_back(std::string(glob_result.gl_pathv[0]));
//        fileNumber ++;
//        searchName = directoryPath + filePrefix + std::to_string(fileNumber) + fileSuffix;
//        glob(searchName.c_str(),GLOB_TILDE,NULL,&glob_result);
////            std::cout << "Number of files " << glob_result.gl_pathc << std::endl;
//    }

    unsigned long pointIndex(0);

    for  (unsigned long i(0); i < pathsToImages.size(); i++) {
        float tmpValue;

        std::vector<float> tmpVec;
        readImageAsVector(pathsToImages[i],tmpVec);

        SquaredEuclideanPoint tmpPoint(tmpVec);
        pointsVec.push_back(tmpPoint);
        pointIndex++;

        if (pointIndex%2000 == 999){
            std::cout << pointIndex+1 << " points read." << std::endl;
        }
    }

    std::cout << "Reading time (ms)" << 1000 * (clock() - timeRead) / CLOCKS_PER_SEC << std::endl;

    //Parallelize
    std::chrono::system_clock::time_point loopTimeStart = std::chrono::system_clock::now();
    std::vector<std::thread> initThreads(numCores);
    unsigned long chunkSize = (maxNumberOfPoints/numCores);

//    #pragma omp parallel for
    for(unsigned t = 0; t < numCores; t++) {

        unsigned long mainPointIndexStart = t * chunkSize;
        unsigned long mainPointIndexEnd = (t + 1) * chunkSize;
        initThreads[t] = std::thread(singleRun, std::ref(pointsVec), t, mainPointIndexStart, mainPointIndexEnd, numberOfInitialPulls, delta);
//        singleRun(pointsVec, mainPointIndexStart, mainPointIndexEnd, numberOfInitialPulls, delta);


    }
    for(unsigned t = 0; t < numCores; t++) {
        initThreads[t].join();
    }
    std::chrono::system_clock::time_point loopTimeEnd = std::chrono::system_clock::now();

    std::cout << "Average time (ms)"
              << std::chrono::duration_cast<std::chrono::milliseconds>(loopTimeEnd - loopTimeStart).count()/maxNumberOfPoints << std::endl;

    return 0;
}