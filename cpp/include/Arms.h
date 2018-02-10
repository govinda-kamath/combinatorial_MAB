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
    float sumOfSquaresOfPulls;
    unsigned long dimension;
    unsigned log10Dimension;
    const templatePoint *point;
    unsigned long id;

    Arm(){
        numberOfPulls = 0;
        sumOfPulls = 0.0;
        upperConfidenceBound = INFINITY;
        lowerConfidenceBound = -INFINITY;
        sumOfSquaresOfPulls = 0.0;
        estimateOfMean = NAN;
        estimateOfSecondMoment = NAN;
    }

    Arm(unsigned long armNumber, const templatePoint &p, unsigned long d) : Arm() {
        id = armNumber;
        point = &p;
        dimension = d;
        log10Dimension = (unsigned) std::ceil(std::log10(dimension));

    }

    ~Arm(){
//        delete point; //no clue why this errors out
    }

    void printArm(){
        std::cout << "Number of pulls" << numberOfPulls*log10Dimension
                  << std::endl;
    }

    friend bool operator> (const Arm& l, const Arm& r)
    {
        return l.lowerConfidenceBound > r.lowerConfidenceBound;
    }

    void updateConfidenceIntervals(float globalSigma, float logDeltaInverse){

        float localSigma, intervalWidth;
        localSigma = globalSigma*100; //Todo: update sigma to new local value
//        localSigma = std::sqrt((sumOfSquaresOfPulls/numberOfPulls -
//                                std::pow(sumOfPulls/numberOfPulls,2)));
        intervalWidth = std::sqrt((localSigma * localSigma * logDeltaInverse)/numberOfPulls);
        upperConfidenceBound = estimateOfMean + intervalWidth;
        lowerConfidenceBound = std::max((float)0.0, estimateOfMean - intervalWidth);
    }

    float pullArm(const templatePoint &p1, float globalSigma,
                  float logDeltaInverse, bool update = true) {
        float sample(-INFINITY);

        if (numberOfPulls >= dimension){
            sample = trueMean();
            numberOfPulls += dimension;
            estimateOfMean = sample;
            upperConfidenceBound = estimateOfMean;
            lowerConfidenceBound = estimateOfMean;
            sumOfPulls += sample*dimension;
            sumOfSquaresOfPulls += std::pow(sample,2)*dimension;
            estimateOfSecondMoment = sample*sample;
        }
        else {
            for(unsigned t = 0; t < log10Dimension*1000000 ; t++)
            {
                sample = point->sampledDistance(p1);
                numberOfPulls++;
                sumOfPulls += sample;
                sumOfSquaresOfPulls += sample * sample;
            }
            estimateOfMean = sumOfPulls / numberOfPulls;
            estimateOfSecondMoment = sumOfSquaresOfPulls / numberOfPulls;
            if (update)
                updateConfidenceIntervals(globalSigma, logDeltaInverse);
        }
        return sample;
    }

    float trueMean(const templatePoint &p1){
        return point->distance(p1)/p1.getVecSize();
    }

    // The child class will has to extend this
    virtual float trueMean(){
        std::cout<<"True mean should not be calculated here!!" <<std::endl;
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
    ArmKNN(unsigned long id, const templatePoint &p) : Arm<templatePoint>(id, p) {}

    ArmKNN(unsigned long id, const templatePoint &p, const templatePoint &fixPoint) : Arm<templatePoint>(id, p, p.vecSize) {
        fixedPoint = &fixPoint;
    }

    using Arm<templatePoint>::pullArm;
    float pullArm(float globalSigma, float logDeltaInverse, bool update = true){
        return pullArm(*fixedPoint, globalSigma, logDeltaInverse, update);
    }

    using Arm<templatePoint>::trueMean;
    float trueMean(){
        return trueMean(*fixedPoint);
    }
};



template <class templatePoint>
class ArmMedoid : public Arm<templatePoint>{
    /*
     * Arms for Medoid*/
public:
    const std::vector<templatePoint> *pointsVec;
    unsigned long numberOfPoints;

    ArmMedoid(unsigned long id, const templatePoint &p) : Arm<templatePoint>(id, p) {}

    ArmMedoid(unsigned long id, const templatePoint &p, const std::vector<templatePoint> &pVec) :
            Arm<templatePoint>(id, p, p.vecSize*pVec.size()) {
        pointsVec = &pVec;
        numberOfPoints = pVec.size();
    }

    using Arm<templatePoint>::pullArm;
    float pullArm(float globalSigma, float logDeltaInverse, bool update = true){
        //Choose a random point
        unsigned long randomCoOrdinate;
        randomCoOrdinate = std::rand() % numberOfPoints;
        return pullArm((*pointsVec)[randomCoOrdinate], globalSigma, logDeltaInverse, update);
    }

    using Arm<templatePoint>::trueMean;
    float trueMean(){
        double Mean = 0;
        for(unsigned long i = 0; i<numberOfPoints; i++){
            Mean += trueMean((*pointsVec)[i]);
        }
        return (float) Mean/((double)numberOfPoints);
    }
};

#endif //COMBINATORIAL_MAB_ARMS_H