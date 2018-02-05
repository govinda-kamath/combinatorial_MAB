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
#include <stdexcept>

#include "utils.h"

/*
 * Finds the nearest k neighbhours of pointsVecL(eft) in pointsVecR(ight)
 */
template <class templatePoint>
class Knn{
public:
    std::vector<templatePoint>& pointsVecL;
    std::vector<templatePoint>& pointsVecR;
    unsigned k;
    unsigned numberOfInitialPulls;
    float delta;

    std::vector<std::vector<ArmKNN<templatePoint>> > nearestNeighbhours;
    std::vector<bool> nnEvaluated;
    std::vector<float> avgNumberOfPulls; //Stat
    bool LEqualR = false; // True when left and right points are the same

    Knn( std::vector<templatePoint>& pVecL, std::vector<templatePoint>& pVecR,
         unsigned NumberOfNeighbhours, unsigned nOfInitialPulls, float deltaAccuracy )
            : pointsVecL(pVecL),
              pointsVecR(pVecR)
    {
        k = NumberOfNeighbhours;
        numberOfInitialPulls = nOfInitialPulls;
        delta = deltaAccuracy;
        nearestNeighbhours.reserve(pVecL.size());
        avgNumberOfPulls.reserve(pVecL.size());
        nnEvaluated.reserve(pVecL.size());
        for(unsigned long i(0); i<pVecL.size(); i++){
            nnEvaluated[i] = false;
            nearestNeighbhours.push_back(std::vector<ArmKNN<templatePoint>>()); //Todo: Bad Code
        }


    }

    Knn( std::vector<templatePoint>& pVecL,
         unsigned NumberOfNeighbhours, unsigned nOfInitialPulls, float deltaAccuracy ) :
            Knn( pVecL, pVecL, NumberOfNeighbhours,  nOfInitialPulls,  deltaAccuracy ) {
        LEqualR = true;
    }


    void run(std::vector<unsigned long> indices){

        for (unsigned long i = 0; i< indices.size(); i++){
            unsigned long index = indices[i];
            if (nnEvaluated[index])
                continue;

            std::vector<ArmKNN<templatePoint> > armsVec;
            for (unsigned long j(0); j < pointsVecR.size(); j++) {
                if ( (j == index) && (LEqualR) ) // Skip the point itself
                    continue;
                ArmKNN<templatePoint> tmpArm(j, pointsVecR[j], pointsVecL[index]);
                armsVec.push_back(tmpArm);
            }

            UCB<ArmKNN<templatePoint> > UCB1(armsVec, delta, k);

            UCB1.initialise(numberOfInitialPulls);

            UCB1.runUCB(200*pointsVecR.size());

            avgNumberOfPulls[index] = UCB1.globalNumberOfPulls/UCB1.numberOfArms;
            nearestNeighbhours[index] = UCB1.topKArms;
            nnEvaluated[index] = true;
        }
    }

    void run(){
        std::vector<unsigned long> indices(pointsVecL.size());
        std::iota(indices.begin(), indices.end(), 0);
        run(indices);
    }

    std::vector<std::vector<ArmKNN<templatePoint>> > get(std::vector<unsigned long> indices){
        std::vector<std::vector<ArmKNN<templatePoint>> > answers;
        for (unsigned long i = 0; i< indices.size(); i++) {
            unsigned long index = indices[i];
            if (!nnEvaluated[index])
                throw std::invalid_argument("Trying to get nearest neighbhour to a point! You have not ran nn on it");
            answers.push_back(nearestNeighbhours[index]);
        }
        return answers;

    }


    void saveAnswers( std::string saveFilePath){

        std::ofstream saveFile;
        saveFile.open (saveFilePath, std::ofstream::out | std::ofstream::app);

        for (unsigned long index = 0; index<pointsVecL.size() ; index++) {
            std::vector<ArmKNN<templatePoint>> topKArms = nearestNeighbhours[index];
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

            saveFile << index << "A";
            for (unsigned i = 0; i < k*5; i++) {
                saveFile <<  " " << topKArmsArgSort[i];
            }
            saveFile << " Av:" << avgNumberOfPulls[index] << "\n";
        }

    }

};