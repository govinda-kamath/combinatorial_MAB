//
// Created by Vivek Kumar Bagaria on 2/3/18.
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
#include <dlib/image_io.h>
#include <dlib/image_transforms.h>
#include "Points.h"
#include "Arms.h"
#include "UCB.h"
#include "UCB_dynamic.h"
#include <stdexcept>
#include "utils.h"

/*
 * Finds the nearest k neighbhours of pointsVectorLeft in pointsVectorRight
 */
template <class templatePoint>
class Knn{
public:
    std::vector<templatePoint> pointsVectorLeft;
    std::vector<templatePoint> pointsVectorRight;
    unsigned k; // Number of nearest neighbhours
    unsigned numberOfInitialPulls; // UCB Parameters
    float delta; // UCB Parameters

    std::vector<std::vector<ArmKNN<templatePoint>> > nearestNeighbours;
    std::vector<short int> nearestNeighboursEvaluated;
    std::vector<float> avgNumberOfPulls; //Stat
    bool leftEqualsRight = false; // True when left and right points are the same

    Knn( const std::vector<templatePoint> &pVecL, const std::vector<templatePoint> &pVecR,
         unsigned NumberOfNeighbours, unsigned noOfInitialPulls, float deltaAccuracy ) {
        pointsVectorLeft = pVecL;
        pointsVectorRight = pVecR;
        initialiseKNN(NumberOfNeighbours, noOfInitialPulls,  deltaAccuracy );
    }

    Knn( const std::vector<templatePoint> &pVecL,
         unsigned NumberOfNeighbours, unsigned noOfInitialPulls, float deltaAccuracy ) {
        pointsVectorLeft = pVecL;
        pointsVectorRight = pVecL;
        leftEqualsRight = true;
        initialiseKNN(NumberOfNeighbours, noOfInitialPulls,  deltaAccuracy );
    }



    void initialiseKNN(unsigned NumberOfNeighbhours, unsigned noOfInitialPulls, float deltaAccuracy ){
        k = NumberOfNeighbhours;
        numberOfInitialPulls = noOfInitialPulls;
        delta = deltaAccuracy;

        nearestNeighbours.reserve(pointsVectorLeft.size());
        avgNumberOfPulls.reserve(pointsVectorLeft.size());

        for(unsigned long i(0); i< pointsVectorLeft.size(); i++){
            nearestNeighboursEvaluated.push_back(false);
            nearestNeighbours.push_back(std::vector<ArmKNN<templatePoint>>()); //Todo: Bad Code
        }

    }

    void run(const std::vector<unsigned long> &indices){

        for (unsigned long i = 0; i < indices.size(); i++){
            unsigned long index = indices[i];
            if (nearestNeighboursEvaluated[index])
                continue;

            std::vector<ArmKNN<templatePoint> > armsVec;
            for (unsigned long j(0); j < pointsVectorRight.size(); j++) {
                if ( (j == index) && (leftEqualsRight) ) // Skip the point itself
                    continue;
                ArmKNN<templatePoint> tmpArm(j, pointsVectorRight[j], pointsVectorLeft[index]);
                armsVec.push_back(tmpArm);
            }

            UCBDynamic<ArmKNN<templatePoint> > UCB1(armsVec, delta, k);

            std::chrono::system_clock::time_point timeStart = std::chrono::system_clock::now();
            UCB1.initialise(numberOfInitialPulls);
            std::chrono::system_clock::time_point timeRunStart = std::chrono::system_clock::now();
            UCB1.runUCB(20000*pointsVectorRight.size());
            std::chrono::system_clock::time_point timeRunEnd = std::chrono::system_clock::now();
//            ArmKNN<templatePoint> tmp = UCB1.bestArm();
//            std::chrono::system_clock::time_point timeTrueMeanEnd = std::chrono::system_clock::now();

            long long int initTime = std::chrono::duration_cast<std::chrono::milliseconds>
                    (timeRunStart - timeStart).count();
            long long int runTime = std::chrono::duration_cast<std::chrono::milliseconds>
                    (timeRunEnd - timeRunStart).count();
//            long long int trueMeanTime = std::chrono::duration_cast<std::chrono::milliseconds>
//                    (timeTrueMeanEnd - timeRunEnd).count();


            avgNumberOfPulls[index] = UCB1.globalNumberOfPulls/UCB1.numberOfArms;
            std::cout << "index " << indices[i] << " Avg Pulls " <<  avgNumberOfPulls[index]
                      << " init time " << initTime << " ms"
                      << " run time " << runTime << " ms"
//                    << " true mean time " << trueMeanTime << " ms"
                      << std::endl;

            nearestNeighbours[index] = UCB1.topKArms;
            nearestNeighboursEvaluated[index] = true;
        }
    }

    void run(){
        std::vector<unsigned long> indices(pointsVectorLeft.size());
        std::iota(indices.begin(), indices.end(), 0);
        run(indices);
    }

    std::vector<std::vector<ArmKNN<templatePoint>> > get(std::vector<unsigned long> indices){
        std::vector<std::vector<ArmKNN<templatePoint>> > answers;
        for (unsigned long i = 0; i< indices.size(); i++) {
            unsigned long index = indices[i];
            if (not nearestNeighboursEvaluated[index])
                throw std::invalid_argument("Trying to get nearest neighbhour to a point! You have not ran nn on it");
            answers.push_back(nearestNeighbours[index]);
        }
        return answers;

    }


    void saveAnswers(std::string saveFilePath){

        std::ofstream saveFile;
        saveFile.open (saveFilePath, std::ofstream::out | std::ofstream::app);

        for (unsigned long index = 0; index<pointsVectorLeft.size() ; index++) {
            if (not nearestNeighboursEvaluated[index])
                continue;
            std::vector<ArmKNN<templatePoint>> topKArms = nearestNeighbours[index];
            std::vector<float> topKArmsTrueMean(k*5);

            saveFile << index << "L ";
            for (unsigned i = 0; i < k*5; i++) {
                saveFile << topKArms[i].id << " ";
                topKArmsTrueMean[i] = topKArms[i].trueMean();
            }
            saveFile << std::endl;

            std::vector<int> topKArmsArgSort(k*5);
            std::iota(topKArmsArgSort.begin(), topKArmsArgSort.end(), 0);
            auto comparator = [&topKArmsTrueMean](int a, int b){ return topKArmsTrueMean[a] < topKArmsTrueMean[b]; };
            std::sort(topKArmsArgSort.begin(), topKArmsArgSort.end(), comparator);

            saveFile << index << "A" << ": ";
            std::cout << "index " << index << ": ";
            for (unsigned i = 0; i < k*5; i++) {
                saveFile <<  " " << topKArmsArgSort[i];
                std::cout <<  " " << topKArmsArgSort[i];
            }
            saveFile << " Av:" << avgNumberOfPulls[index] << "\n";
            std::cout << " Av:" << avgNumberOfPulls[index] << "\n";


            std::cout << "index " << index << ": " ;
            for (unsigned i = 0; i < k*5; i++) {
                saveFile <<  " " << topKArmsArgSort[i];
                std::cout <<  " " << topKArms[i].id;
            }
            saveFile << " Av:" << avgNumberOfPulls[index] << "\n";
            std::cout << " Av:" << avgNumberOfPulls[index] << "\n";
        }
    }
};