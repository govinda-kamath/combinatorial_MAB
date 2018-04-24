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
    unsigned int sampleSize;

    //Outputs
    std::vector<std::vector<ArmKNN<templatePoint>> > nearestNeighbours;
    std::vector< std::vector<unsigned long> > finalNumberOfPulls;
    std::vector< std::vector<unsigned long> > finalSortedOrder;
    std::vector<std::vector<ArmKNN<templatePoint>> > nearestNeighboursBrute;
    std::vector<short int> nearestNeighboursEvaluated;
    std::vector<long long int> initTime;
    std::vector<long long int> runTime;

    std::vector<float> avgNumberOfPulls; //Statistics
    bool leftEqualsRight = false; // True when left and right points are the same

    Knn( const std::vector<templatePoint> &pVecL, const std::vector<templatePoint> &pVecR,
         unsigned NumberOfNeighbours, unsigned noOfInitialPulls, float deltaAccuracy, unsigned int sSize ) {
        pointsVectorLeft = pVecL;
        pointsVectorRight = pVecR;
        sampleSize = sSize;
        initialiseKNN(NumberOfNeighbours, noOfInitialPulls,  deltaAccuracy );
    }

    Knn( const std::vector<templatePoint> &pVecL,
         unsigned NumberOfNeighbours, unsigned noOfInitialPulls, float deltaAccuracy , unsigned int sSize ) {
        pointsVectorLeft = pVecL;
        pointsVectorRight = pVecL;
        leftEqualsRight = true;
        sampleSize = sSize;
        initialiseKNN(NumberOfNeighbours, noOfInitialPulls,  deltaAccuracy );
    }



    void initialiseKNN(unsigned NumberOfNeighbhours, unsigned noOfInitialPulls, float deltaAccuracy ){
        k = NumberOfNeighbhours;
        numberOfInitialPulls = noOfInitialPulls;
        delta = deltaAccuracy;

        nearestNeighbours.reserve(pointsVectorLeft.size());
        nearestNeighboursBrute.reserve(pointsVectorLeft.size());
        avgNumberOfPulls.reserve(pointsVectorLeft.size());

        for(unsigned long i(0); i< pointsVectorLeft.size(); i++){
            nearestNeighboursEvaluated.push_back(false);
            nearestNeighbours.push_back(std::vector<ArmKNN<templatePoint>>()); //Todo: Bad Code
            finalNumberOfPulls.push_back(std::vector<unsigned long>()); //Todo: Bad Code
            finalSortedOrder.push_back(std::vector<unsigned long>()); //Todo: Bad Code
#ifdef Brute
            nearestNeighboursBrute.push_back(std::vector<ArmKNN<templatePoint>>()); //Todo: Bad Code
#endif
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

            UCBDynamic<ArmKNN<templatePoint> > UCB1(armsVec, delta, k, 5*k, sampleSize);
            std::chrono::system_clock::time_point timeStart = std::chrono::system_clock::now();
//#define Brute
//#ifndef Brute
            UCB1.initialise(numberOfInitialPulls);
//            std::chrono::system_clock::time_point timeRunStart = std::chrono::system_clock::now();
            UCB1.runUCB(2000*pointsVectorRight.size());
//#endif

// Stats
#ifdef Brute
            std::cout << "Running Brute" << std::endl;
//            std::chrono::system_clock::time_point timeRunStart = std::chrono::system_clock::now();
            UCB1.armsKeepFromArmsContainerBrute();
#endif
            std::chrono::system_clock::time_point timeRunEnd = std::chrono::system_clock::now();

            initTime.push_back(std::chrono::duration_cast<std::chrono::milliseconds>
                    (timeRunStart - timeStart).count());
            runTime.push_back(std::chrono::duration_cast<std::chrono::milliseconds>
                    (timeRunEnd - timeRunStart).count());


            UCB1.storeExtraTopArms();
            avgNumberOfPulls[index] = UCB1.globalNumberOfPulls/UCB1.numberOfArms;
//            if (index%25==0){
                std::cout << "index " << indices[i] << " Avg Pulls " <<  avgNumberOfPulls[index]
                          << " init time " << initTime[index] << " ms"
                          << " run time " << runTime[index] << " ms"
                          << std::endl;
//            }
            nearestNeighbours[index] = UCB1.topKArms;
            finalSortedOrder[index] = UCB1.finalSortedOrder;
            finalNumberOfPulls[index] = UCB1.finalNumberOfPulls;
#ifdef Brute
            nearestNeighboursBrute[index] = UCB1.bruteBestArms();
#endif

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
                throw std::invalid_argument("Trying to get nearest neighbhour to a point! You havent yet run NN on it");
            answers.push_back(nearestNeighbours[index]);
        }
        return answers;

    }


    void saveAnswers(std::string saveFolderPath){

        std::string  n = std::to_string(pointsVectorLeft.size());
        std::string  d = std::to_string(nearestNeighbours[0][0].dimension);
        int sum = 0 ;
        int total = 0;
        for (unsigned long index = 0; index<pointsVectorLeft.size() ; index++) {
            if (not nearestNeighboursEvaluated[index]){
                total += 1;
                continue;
            }

            std::string saveFilePath = saveFolderPath+"n_"+n+"_d_"+d+"_k_"+std::to_string(k)+"_index_"+std::to_string(index);
//            std::cout << saveFilePath << std::endl;
            std::ofstream saveFile;
            saveFile.open (saveFilePath, std::ofstream::out | std::ofstream::app);

            saveFile << "AveragePulls," << avgNumberOfPulls[index] << "\n";
            saveFile << "InitTime," << initTime[index] << "\n";
            saveFile << "RunTime," << runTime[index] << "\n";
            saveFile << "NumberOfInitialPulls," << numberOfInitialPulls << "\n";
            saveFile << "SampleSize," << sampleSize << "\n";
            saveFile << "n," << n << "\n";
            saveFile << "d," << d << "\n";
            saveFile << "k," << k << "\n";

            std::vector<ArmKNN<templatePoint>> topKArms = nearestNeighbours[index];
#ifdef Brute
            std::vector<ArmKNN<templatePoint>> topKArmsBrute = nearestNeighboursBrute[index];
#endif
            std::vector<float> topKArmsTrueMean(k*5);
            std::vector<float> topKArmsTrueMeanBrute(k*5);


            saveFile << "Answer,";
            for (unsigned i = 0; i < k*5; i++) {
                saveFile << topKArms[i].id << ",";
#ifdef Brute
                topKArmsTrueMeanBrute[i] = topKArmsBrute[i].trueMean();
#endif
                topKArmsTrueMean[i] = topKArms[i].trueMean();
            }
#ifndef Brute

            saveFile << std::endl;

            std::vector<int> topKArmsArgSort(k*5);
            std::iota(topKArmsArgSort.begin(), topKArmsArgSort.end(), 0);
            auto comparator = [&topKArmsTrueMean](int a, int b){ return topKArmsTrueMean[a] < topKArmsTrueMean[b]; };
            std::sort(topKArmsArgSort.begin(), topKArmsArgSort.end(), comparator);
            std::cout << "\nPoint number " << index << ": Rank:\t";

            saveFile << "Position";
            for (unsigned i = 0; i < k*5; i++) {
                saveFile <<  "," << topKArmsArgSort[i];
                std::cout <<  " " << topKArmsArgSort[i];
            }
            saveFile << std::endl;

            saveFile << "AllPullsNumber";
            for (unsigned i = 0; i < pointsVectorRight.size(); i++) {
                saveFile <<  "," << finalNumberOfPulls[index][i];
            }
            saveFile << std::endl;

            saveFile << "AllPullsIndex";
            for (unsigned i = 0; i < pointsVectorRight.size(); i++) {
                saveFile <<  "," << finalSortedOrder[index][i];
            }
            saveFile << std::endl;



            std::cout << " Av:" << avgNumberOfPulls[index] << "\n";


            std::cout << "\nPoint number " << index << ": ArmId:\t";
            for (unsigned i = 0; i < k*5; i++) {
//                saveFile <<  " " << topKArmsArgSort[i];
                std::cout <<  " " << topKArms[i].id;
            }

#endif
#ifdef Brute
            bool flag = true;
            for (unsigned i = 0; i < k; i++) {
                if (topKArms[i].id!=topKArmsBrute[i].id)
                    flag = false;
            }
//
//
////            saveFile << " Av:" << avgNumberOfPulls[index] << "\n";
//            std::cout << " Av:" << avgNumberOfPulls[index];
            std::cout << "index " << index << ": ";

            std::cout << "Verdict: " << flag << "\n";
            if (flag == false){
                std::cout << "\nUCB: " ;
                for (unsigned i = 0; i < 2*k; i++) {
                    std::cout << topKArms[i].id << " ";
                }

                std::cout << "\nBrute: " ;
                for (unsigned i = 0; i < 2*k; i++) {
                    std::cout << topKArmsBrute[i].id << " ";
                }

                for (unsigned i = 0; i < 2*k; i++) {
                    std::cout << topKArmsTrueMean[i] << " ";
                }
                std::cout << "\nBrute: " ;
                for (unsigned i = 0; i < 2*k; i++) {
                    std::cout << topKArmsTrueMeanBrute[i] << " ";
                }
            }
            sum += (int) flag;
#endif
        }

#ifdef Brute
        std::cout << "Total " << pointsVectorLeft.size() - total
                              << "Correct " << sum << std::endl;
#endif
    }
};