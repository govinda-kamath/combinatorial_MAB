//
// Created by Govinda Kamath on 2/5/18.
//

#ifndef COMBINATORIAL_MAB_KMEANS_H
#define COMBINATORIAL_MAB_KMEANS_H

#endif //COMBINATORIAL_MAB_KMEANS_H
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
#include "../deprecated/UCB.h"
#include <stdexcept>
#include "utils.h"
#include "Knn.h"

template <class templatePoint>
class Kmeans {
public:

    unsigned numberOfInitialPulls; // UCB Parameters
    float delta; //  UCB Parameters
    unsigned int sampleSize;
    unsigned int k;
    unsigned long n;
    std::string savefilePath;

    std::vector< std::vector<unsigned long> > clusters;
    std::vector<templatePoint> centersVec;
    std::vector<templatePoint> pointsVec;

    //Output
    std::vector<float> finalNumberOfPulls;
    float totalTime;
    float avgNumPulls;

    Kmeans( std::vector<templatePoint>& pVec, std::vector<templatePoint>& cVec,
            unsigned noOfInitialPulls, float deltaAccuracy , unsigned int sSize,
            std::string sFilePath){
        centersVec = cVec;
        pointsVec = pVec;
        numberOfInitialPulls = noOfInitialPulls;
        delta = deltaAccuracy;
        sampleSize = sSize;
        k = cVec.size();
        n = pVec.size();
        savefilePath = sFilePath;
    }



    void maximization(){
        for(unsigned j(0); j<k ;j++){
            clusters.push_back(std::vector<unsigned long>());
        }
        Knn<templatePoint> knn( pointsVec, centersVec, 1, numberOfInitialPulls, delta, sampleSize, "/tmp");
        std::chrono::system_clock::time_point timeStart = std::chrono::system_clock::now();
        knn.run();
        std::chrono::system_clock::time_point timeEnd = std::chrono::system_clock::now();

        for (unsigned long i(0); i < pointsVec.size(); i ++){
            unsigned long index = knn.allNearestNeighbhourIds[i][0];
            clusters[index].push_back(i);
        }
        finalNumberOfPulls = knn.avgNumberOfPulls;
        totalTime = std::chrono::duration_cast<std::chrono::microseconds>
                (timeStart - timeEnd).count();



    }

    // Todo: To test the function
    void expectation(){
        for(unsigned i(0); i < k ; i++){
            std::vector<float> newCenter(n, 0.0);
            for(unsigned long j(0); j<clusters[i].size(); j++){
                std::transform(newCenter.begin(), newCenter.end(), pointsVec[clusters[i][j]].point.begin(),
                               std::back_inserter(newCenter), std::plus<float>());
            }

            std::transform( newCenter.begin(), newCenter.end(), newCenter.begin(),
                            std::bind2nd(std::divides<float>(), clusters[i].size()) );
            centersVec[i].point = newCenter;
        }
    }

    void saveAnswers() {
        //Variables
        std::string d = std::to_string(pointsVec[0].vecSize);
        std::string saveFilePath =
                savefilePath + "n_" + std::to_string(n) + "_d_" + d + "_k_" + std::to_string(k);
        std::ofstream saveFile;
        saveFile.open(saveFilePath, std::ofstream::out | std::ofstream::trunc);
        std::cout << saveFilePath << " PATH " << finalNumberOfPulls.size() << std::endl;
        avgNumPulls = 0;
        for (unsigned long i(0); i < n; i++){
            avgNumPulls += finalNumberOfPulls[i];
        }
        avgNumPulls =  avgNumPulls/n;

        //Save single stats
        saveFile << "AveragePulls," << avgNumPulls << std::endl;
        saveFile << "RunTime," << totalTime << std::endl;
        saveFile << "NumberOfInitialPulls," << numberOfInitialPulls<< std::endl;
        saveFile << "SampleSize," << sampleSize << std::endl;
        saveFile << "n," << n << std::endl;
        saveFile << "d," << d << std::endl;
        saveFile << "k," << k << std::endl;
        saveFile << "AllPullsNumber";

        for (unsigned i(0); i < n; i++) {
            saveFile <<  "," << finalNumberOfPulls[i];
        }
        saveFile << std::endl;

    }



};
