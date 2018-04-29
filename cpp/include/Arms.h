//
// Created by Vivek Kumar Bagaria on 1/31/18.
//

#ifndef COMBINATORIAL_MAB_ARMS_H
#define COMBINATORIAL_MAB_ARMS_H

#include "Points.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <thread>

template <class templatePoint>
class Arm {
    /* A base class of arms */
public:
    unsigned long numberOfPulls;
    float sumOfPulls;
    float upperConfidenceBound;
    float lowerConfidenceBound;
    float estimateOfMean;
    float estimateOfSecondMoment;
    float localSigma;
    float sumOfSquaresOfPulls;
    unsigned long dimension;
    unsigned log10Dimension;
    std::unordered_map<std::string, float> misc;
    const templatePoint *point;
    unsigned long id;
    float trueMeanValue;

    Arm(){
        numberOfPulls = 0;
        sumOfPulls = 0.0;
        upperConfidenceBound = INFINITY;
        lowerConfidenceBound = -INFINITY;
        sumOfSquaresOfPulls = 0.0;
        estimateOfMean = NAN;
        estimateOfSecondMoment = NAN;
        localSigma = NAN;
        trueMeanValue = INFINITY;
    }

    Arm(unsigned long armNumber, unsigned long d) : Arm() {
        id = armNumber;
        dimension = d;
        log10Dimension = (unsigned) std::ceil(std::log10(dimension));
    }

    Arm(unsigned long armNumber, const templatePoint &p, unsigned long d) : Arm() {
        id = armNumber;
        point = &p;
        dimension = d;
        log10Dimension = (unsigned) std::ceil(std::log10(dimension));
//        std::cout << "log10D " << log10Dimension << std::endl;

    }

    ~Arm(){
//        delete point; //no clue why this errors out
    }

    friend bool operator> (const Arm &l, const Arm &r)
    {
        return l.lowerConfidenceBound > r.lowerConfidenceBound;
    }

    void updateConfidenceIntervals(float globalSigma, unsigned long long globalNumberOfPulls, float logDeltaInverse){


        float compositeSigma, intervalWidth;
//        compositeSigma = globalSigma; //Todo: update sigma to new local value
        localSigma = std::sqrt(estimateOfSecondMoment - std::pow(estimateOfMean,2));
        if (localSigma<0){
            std::cout << "duck!" <<std::endl;
        }
        float frac = numberOfPulls/dimension;
        if (frac>=1){
            frac = 1;
        }
//        frac = 0;
        compositeSigma = std::sqrt( localSigma*localSigma*frac +  globalSigma*globalSigma*(1- frac));


        intervalWidth = std::sqrt((compositeSigma * compositeSigma * logDeltaInverse)/(float)numberOfPulls);
        upperConfidenceBound = estimateOfMean + intervalWidth;
        lowerConfidenceBound = std::max((float)0.0, estimateOfMean - intervalWidth);
    }


    void updateMeanAndSecondMoment(){
        if (trueMeanValue != INFINITY){
            estimateOfMean = trueMeanValue;
            estimateOfSecondMoment = trueMeanValue*trueMeanValue;
        }
        else{
            estimateOfMean = sumOfPulls/numberOfPulls;
            estimateOfSecondMoment = sumOfSquaresOfPulls/numberOfPulls;
        }
    }

    void warmInitialise(unsigned long numArmPulls, float armSumOfPulls, float armSumOfSquaresOfPulls,
                        float tMeanValue){
        numberOfPulls = numArmPulls;
        sumOfPulls = armSumOfPulls;
        sumOfSquaresOfPulls = armSumOfSquaresOfPulls;
        trueMeanValue = tMeanValue;
        updateMeanAndSecondMoment();
    }

    virtual std::pair<float, float> pullArm(const templatePoint &p1, const templatePoint &p2, float globalSigma,
                  unsigned long long globalNumberOfPulls,  float logDeltaInverse, bool update, unsigned sampleSize, float LCBofSecondBestArm){
        std::pair<float, float> sample;
        if (numberOfPulls >= dimension){
            float tMean(-INFINITY);
            tMean = trueMean();
            numberOfPulls += sampleSize;
            estimateOfMean = tMean;
            upperConfidenceBound = estimateOfMean;
            lowerConfidenceBound = estimateOfMean;
            sumOfPulls += tMean*dimension;
            sumOfSquaresOfPulls += std::pow(tMean,2)*dimension;
            estimateOfSecondMoment = tMean*tMean;
            sample.first = tMean;
            sample.second = tMean*tMean;
        }
        else {
//            unsigned tmp = 0;
           do {
//               tmp++;
//               std::cout << id << "\t" << "tmp = " << tmp << "\t"
//                         << lowerConfidenceBound << "\t"
//                         << "\t" << LCBofSecondBestArm
//                         << "\t" << upperConfidenceBound
//                         << std::endl;
               sample = p2.sampleDistance(p1, sampleSize);
               numberOfPulls += sampleSize;
               sumOfPulls += sample.first;
               sumOfSquaresOfPulls += sample.second;
               estimateOfMean = sumOfPulls / numberOfPulls;
               estimateOfSecondMoment = sumOfSquaresOfPulls / numberOfPulls;
               if (update)
                   updateConfidenceIntervals(globalSigma, globalNumberOfPulls, logDeltaInverse);
           } while (lowerConfidenceBound < LCBofSecondBestArm && upperConfidenceBound  > LCBofSecondBestArm && update && false);
        }
        return sample;
    }

    virtual std::pair<float, float> pullArm(const templatePoint &p1, float globalSigma, unsigned long long globalNumberOfPulls,
                  float logDeltaInverse, bool update, unsigned sampleSize, float LCBofSecondBestArm) {
        return pullArm(p1, *point,  globalSigma,globalNumberOfPulls, logDeltaInverse, update, sampleSize, LCBofSecondBestArm);

    }

    float trueMean(const templatePoint &p1){
        if (trueMeanValue == INFINITY){
            trueMeanValue = point->distance(p1)/p1.getVecSize();
        }
        return trueMeanValue;
    }

    // The child class will has to extend this
    virtual float trueMean(){
        std::cout << "True mean should not be calculated here!!" << std::endl;
        return -1;
    };
};


template <class templatePoint>
class ArmKNN : public Arm<templatePoint>{
    /*
     * Arms for k-Nearest Neighbours*/
public:
    const templatePoint *fixedPoint;

    ArmKNN(): Arm<templatePoint>(){}

    using Arm<templatePoint>::id;
    ArmKNN(unsigned long armNumber){
        id = armNumber;
    }

    ArmKNN(unsigned long id, const templatePoint &p) : Arm<templatePoint>(id, p) {}

    ArmKNN(unsigned long id, const templatePoint &p, const templatePoint &fixPoint) : Arm<templatePoint>(id, p, p.vecSize) {
        fixedPoint = &fixPoint;
    }

    using Arm<templatePoint>::pullArm;
    virtual std::pair<float, float> pullArm(float globalSigma, unsigned long long globalNumberOfPulls,
                  float logDeltaInverse, bool update , unsigned sampleSize, float LCBofSecondBestArm){
        return pullArm(*fixedPoint, globalSigma, globalNumberOfPulls, logDeltaInverse, update, sampleSize, LCBofSecondBestArm);
    }

    using Arm<templatePoint>::trueMean;
    float trueMean(){
        return trueMean(*fixedPoint);
    }

    using Arm<templatePoint>::estimateOfMean;
    using Arm<templatePoint>::upperConfidenceBound;
    using Arm<templatePoint>::lowerConfidenceBound;

    std::unordered_map<std::string, float> trueMeanUpdate(){
        float localSumOfPulls = trueMean(*fixedPoint);
        float localSumOfSquaresOfPulls = localSumOfPulls*localSumOfPulls;

        std::unordered_map<std::string, float> result;
        unsigned  long d = 4000; //ToDo: Change this

        result.insert( std::make_pair<std::string, float>("sumOfPulls", localSumOfPulls*d));
        result.insert( std::make_pair<std::string, float>("sumOfSquaresPulls", localSumOfSquaresOfPulls*d));
        result.insert( std::make_pair<std::string, float>("effectiveDimension", (float) d));
        estimateOfMean = localSumOfPulls;
        upperConfidenceBound = estimateOfMean;
        lowerConfidenceBound = estimateOfMean;
        return result;
    }
};



template <class templatePoint>
class ArmMedoid : public Arm<templatePoint>{
    /*
     * Arms for Medoid*/
public:
    const std::vector<templatePoint> *pointsVec;
    unsigned long numberOfPoints;

//    ArmMedoid(unsigned long id, const templatePoint &p) : Arm<templatePoint>(id, p) {}

    ArmMedoid(unsigned long id, const templatePoint &p, const std::vector<templatePoint> &pVec) :
            Arm<templatePoint>(id, p, p.vecSize*pVec.size()) {
        pointsVec = &pVec;
        numberOfPoints = pVec.size();
    }

    using Arm<templatePoint>::pullArm;
    virtual std::pair<float, float> pullArm(float globalSigma, unsigned long long globalNumberOfPulls,
                  float logDeltaInverse, bool update, unsigned sampleSize, float LCBofSecondBestArm){
        //Choose a random point
        unsigned long randomCoOrdinate;
        randomCoOrdinate = std::rand() % numberOfPoints;
        return pullArm((*pointsVec)[randomCoOrdinate], globalSigma, globalNumberOfPulls, logDeltaInverse, update, sampleSize, LCBofSecondBestArm);
    }

    using Arm<templatePoint>::trueMean;
    float trueMean(){
        float mean = 0;
        for(unsigned long i = 0; i<numberOfPoints; i++){
            float tmp = trueMean((*pointsVec)[i]);
            mean += tmp;
        }
        return mean/((float)numberOfPoints);
    }

    using Arm<templatePoint>::estimateOfMean;
    using Arm<templatePoint>::upperConfidenceBound;
    using Arm<templatePoint>::lowerConfidenceBound;
    std::unordered_map<std::string, float> trueMeanUpdate(){
        float localSumOfPulls = 0;
        float localSumOfSquaresOfPulls = 0;
        for(unsigned long i = 0; i<numberOfPoints; i++){
            float colMean = trueMean((*pointsVec)[i]);
            localSumOfPulls += colMean;
            localSumOfSquaresOfPulls += colMean*colMean;
        }
        std::unordered_map<std::string, float> result;
        unsigned  long d = 4000;
        result.insert( std::make_pair<std::string, float>("sumOfPulls", localSumOfPulls*d));
        result.insert( std::make_pair<std::string, float>("sumOfSquaresPulls", localSumOfSquaresOfPulls*d));
        result.insert( std::make_pair<std::string, float>("effectiveDimension", (float) d*numberOfPoints));
        estimateOfMean = localSumOfPulls;
        upperConfidenceBound = estimateOfMean;
        lowerConfidenceBound = estimateOfMean;
        return result;
    }
};

template <class templatePoint>
class ArmHeirarchical : public Arm<GroupPoint<templatePoint> >{
/*
     * Arms for Heirarchical Clustering */
public:
    GroupPoint<templatePoint> *leftGroupPoint, *rightGroupPoint;
    unsigned long leftGroupID, rightGroupID;

    using Arm<GroupPoint<templatePoint> > ::id;
    ArmHeirarchical(unsigned long armNumber){
        id = armNumber;
    }

    ArmHeirarchical(unsigned long id, GroupPoint<templatePoint> &p1, GroupPoint<templatePoint> &p2) :
            Arm<GroupPoint<templatePoint> >(id, p1.getVecSize()){
        leftGroupPoint = &p1;
        rightGroupPoint = &p2;
        leftGroupID = p1.groupID;
        rightGroupID = p2.groupID;
    }

    using Arm<GroupPoint<templatePoint> >::pullArm;
    virtual std::pair<float, float> pullArm(float globalSigma, unsigned long long globalNumberOfPulls,
                  float logDeltaInverse, bool update, unsigned sampleSize, float LCBofSecondBestArm){
        return pullArm(*leftGroupPoint, *rightGroupPoint, globalSigma, globalNumberOfPulls, logDeltaInverse, update, sampleSize, LCBofSecondBestArm);
    }

    using Arm<GroupPoint<templatePoint> >::trueMeanValue;
    using Arm<GroupPoint<templatePoint> >::dimension;
    float trueMean(){
        if (trueMeanValue == INFINITY){
            trueMeanValue = leftGroupPoint->distance(*rightGroupPoint)/(dimension);
        }
//        std::cout<< leftGroupID << " " << rightGroupID << " " << trueMeanValue << std::endl;
        return trueMeanValue;
    }
};

#endif //COMBINATORIAL_MAB_ARMS_H