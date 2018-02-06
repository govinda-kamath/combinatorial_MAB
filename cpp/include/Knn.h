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
 * Finds the nearest k neighbhours of pointsVectorLeft in pointsVectorRight
 */
template <class templatePoint>
class Knn{
public:
    std::vector<templatePoint> pointsVectorLeft;
    std::vector<templatePoint> pointsVectorRight;
    unsigned k;
    unsigned numberOfInitialPulls;
    float delta;

    std::vector<std::vector<ArmKNN<templatePoint>> > nearestNeighbours;
    std::vector<short int> nearestNeighboursEvaluated;
    std::vector<float> avgNumberOfPulls; //Stat
    bool leftEqualsRight = false; // True when left and right points are the same


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


    Knn( std::vector<templatePoint>& pVecL, std::vector<templatePoint>& pVecR,
         unsigned NumberOfNeighbours, unsigned noOfInitialPulls, float deltaAccuracy ) {
        pointsVectorLeft = pVecL;
        pointsVectorRight = pVecR;
        initialiseKNN(NumberOfNeighbours, noOfInitialPulls,  deltaAccuracy );
    }

    Knn( std::vector<templatePoint>& pVecL,
         unsigned NumberOfNeighbours, unsigned noOfInitialPulls, float deltaAccuracy ) {
        pointsVectorLeft = pVecL;
        pointsVectorRight = pVecL;
        leftEqualsRight = true;
        initialiseKNN(NumberOfNeighbours, noOfInitialPulls,  deltaAccuracy );
    }



    void run(std::vector<unsigned long> indices){

        for (int index(0); index < 100; index++) {
            std::cout << nearestNeighboursEvaluated[index] << " ";
        }
        std::cout << std::endl << ":done" << std::endl;

        for (int index(0); index < 100; index++)
            std::cout << nearestNeighboursEvaluated[index] << " ";
        std::cout << std::endl << ":done" << std::endl;


        for (unsigned long i = 0; i < indices.size(); i++){
            std::cout << i << " "<< indices[i] << " " << indices.size() << std::endl;
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

            std::cout << indices[i] << "adf"<< std::endl;
            UCB<ArmKNN<templatePoint> > UCB1(armsVec, delta, k);

            UCB1.initialise(numberOfInitialPulls);

            std::cout << indices[i] << "adfw" << std::endl;
            UCB1.runUCB(200*pointsVectorRight.size());
            std::cout << indices[i] << "adfwdlk" << std::endl;


//            avgNumberOfPulls.push_back(UCB1.globalNumberOfPulls/UCB1.numberOfArms);
//            nearestNeighbours.push_back(UCB1.topKArms);
//            nearestNeighboursEvaluated[index] = true;
            for (int index1(0); index1 < 20; index1++)
                std::cout << nearestNeighboursEvaluated[index1] << " ";
            std::cout << std::endl << ":done" << std::endl;
            std::cout << indices[i] << "adfwdlkadfas" << std::endl;
        }
//        std::cout  << "adfwdlkadfasasd" << std::endl;
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


    void saveAnswers( std::string saveFilePath){

        std::ofstream saveFile;
        saveFile.open (saveFilePath, std::ofstream::out | std::ofstream::app);

        for (unsigned long index = 0; index<pointsVectorLeft.size() ; index++) {
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

            saveFile << index << "A";
            for (unsigned i = 0; i < k*5; i++) {
                saveFile <<  " " << topKArmsArgSort[i];
            }
            saveFile << " Av:" << avgNumberOfPulls[index] << "\n";
        }
    }
};


//template <class templatePoint>
//class Kmeans : public Knn<templatePoint>{
//public:
//Kmeans( std::vector<templatePoint>& pVecL, std::vector<templatePoint>& pVecR,
//unsigned noOfInitialPulls, float deltaAccuracy ):
//Knn<templatePoint>(pVecL,  pVecR, 1,  noOfInitialPulls,  deltaAccuracy ){}
//
//using Knn<templatePoint>::run;
//
//using Knn<templatePoint>::saveAnswers;
//
//};