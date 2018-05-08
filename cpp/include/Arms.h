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
#include <random>
#include <algorithm>

//#include <nanoflann.hpp>

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
    std::unordered_map<std::string, float> misc;
    const templatePoint *point;
    long id;
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
    }

    Arm(unsigned long armNumber, const templatePoint &p, unsigned long d) : Arm() {
        id = armNumber;
        point = &p;
        dimension = d;
    }

    ~Arm(){
//        delete point; //no clue why this errors out
    }

    friend bool operator> (const Arm &l, const Arm &r)
    {
        return l.lowerConfidenceBound > r.lowerConfidenceBound;
    }

    void updateConfidenceIntervals(float globalSigma, unsigned long globalNumberOfPulls, float logDeltaInverse){


        float compositeSigma, intervalWidth;
        localSigma = std::sqrt(estimateOfSecondMoment - std::pow(estimateOfMean,2));
        if (localSigma<0){
            std::cout << "Abort mission!! Fundamental error" <<std::endl;
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
                  unsigned long globalNumberOfPulls,  float logDeltaInverse, bool update, unsigned sampleSize, float LCBofSecondBestArm){
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

    virtual std::pair<float, float> pullArm(const templatePoint &p1, float globalSigma, unsigned long globalNumberOfPulls,
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
    virtual std::pair<float, float> pullArm(float globalSigma, unsigned long globalNumberOfPulls,
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
    virtual std::pair<float, float> pullArm(float globalSigma, unsigned long globalNumberOfPulls,
                  float logDeltaInverse, bool update, unsigned sampleSize, float LCBofSecondBestArm){
        //Choose a random point
        unsigned long randomCoOrdinate;
        randomCoOrdinate = std::rand() % numberOfPoints;
        return pullArm((*pointsVec)[randomCoOrdinate], globalSigma, globalNumberOfPulls, logDeltaInverse,
                       update, sampleSize, LCBofSecondBestArm);
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
    virtual std::pair<float, float> pullArm(float globalSigma, unsigned long globalNumberOfPulls,
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

//H(X_1,...,X_d)
template <class templatePoint>
class ArmEntropyContinuous: public Arm<templatePoint>{
//    std::vector<templatePoint> allPoints;
public:
    std::vector<templatePoint> sampledPoints;
    std::vector<float> nearestNeighbhourDistance;
    std::vector<unsigned long> nearestNeighbhourIndex;
    std::vector<unsigned long> *shuffledRows;
    unsigned long maxSize;

    using Arm<templatePoint>::numberOfPulls;
    using Arm<templatePoint>::sumOfPulls;
    using Arm<templatePoint>::sumOfSquaresOfPulls;
    using Arm<templatePoint>::estimateOfMean;
    using Arm<templatePoint>::estimateOfSecondMoment;
    using Arm<templatePoint>::upperConfidenceBound;
    using Arm<templatePoint>::lowerConfidenceBound;
    using Arm<templatePoint>::localSigma;
    using Arm<templatePoint>::trueMeanValue;

    const float chi_d = 3.0; //the constant multiplying the variance


    ArmEntropyContinuous(unsigned long id, unsigned long maxSize_){
        maxSize = maxSize_;
        for(unsigned long i(0); i< maxSize; i++) {
            nearestNeighbhourIndex.push_back(-1);
            nearestNeighbhourDistance.push_back(INFINITY);
        }
        numberOfPulls = 1;
        trueMeanValue = INFINITY;
    }

    using Arm<templatePoint>::pullArm;
    virtual std::pair<float, float> pullArm(std::vector<templatePoint> samplePoints, float globalSigma,
                                            float logDeltaInverse, bool update, float LCBofSecondBestArm) {
        float sampleSum(0), sampleSquareSum(0);
        for (unsigned i(0); i < samplePoints.size(); i++) {

            templatePoint p = samplePoints[i];
            float minDistance(INFINITY);
            unsigned long minDistIndex(-1);
            //Calculating the nearest neighbhour for the new sampled point.
            for (unsigned j(0); j < numberOfPulls - 1 + i; j++) {
                float dist = p.distance(sampledPoints[j]);
                if (dist < minDistance) {
                    dist = minDistance;
                    minDistIndex = j;
                }

                if (dist < nearestNeighbhourDistance[j]) {
                    sampleSum -= std::log(nearestNeighbhourDistance[j]);
                    sampleSquareSum -= std::pow(std::log(nearestNeighbhourDistance[j]), 2);
                    sampleSum += std::log(dist);
                    sampleSquareSum += std::pow(std::log(dist), 2);

                    nearestNeighbhourDistance[j] = dist;
                    nearestNeighbhourIndex[j] = numberOfPulls + i;
                }
            }
            if (minDistance != INFINITY) {
                sampleSum += std::log(minDistance);
                sampleSquareSum += std::pow(std::log(minDistance), 2);
                nearestNeighbhourDistance[numberOfPulls + i] = minDistance;
                nearestNeighbhourIndex[numberOfPulls + i] = minDistIndex;
                sampledPoints.push_back(p);
            }
        }

        if (sampleSum != INFINITY)
        {
            numberOfPulls += samplePoints.size();
            sumOfPulls += sampleSum;
            sumOfSquaresOfPulls += sampleSquareSum;
            estimateOfMean = sumOfPulls / numberOfPulls;
            estimateOfSecondMoment = sumOfSquaresOfPulls / numberOfPulls;
        }

        if (update){
            updateConfidenceIntervals(globalSigma, logDeltaInverse);
        }
        return std::make_pair(sampleSum, sampleSquareSum);
    }


    using Arm<templatePoint>::updateConfidenceIntervals;
    void updateConfidenceIntervals(float globalSigma, float logDeltaInverse){
        float compositeSigma, intervalWidth;
        localSigma = std::sqrt(estimateOfSecondMoment - std::pow(estimateOfMean,2));
        if (localSigma<0){
            std::cout << "Abort mission!! Fundamental error" <<std::endl;
        }
        float frac = numberOfPulls/maxSize;
        if (frac>=1){
            intervalWidth = 0;
        }
        else{
            compositeSigma = std::sqrt(localSigma*localSigma*frac +  globalSigma*globalSigma*(1- frac)+ chi_d);
            intervalWidth = std::sqrt((compositeSigma * compositeSigma * logDeltaInverse)/(float)numberOfPulls);
        }
        upperConfidenceBound = estimateOfMean + intervalWidth;
        lowerConfidenceBound = std::max((float)0.0, estimateOfMean - intervalWidth);
    }
/*
    // To test
    using Arm<templatePoint>::trueMean;
    float trueMean(const templatePoint &p1){
        if (trueMeanValue == INFINITY){
            float sampleSum(0);
            for(unsigned i(0); i<allPoints.size(); i++){
                float minDistance = INFINITY;
                unsigned long minDistIndex = -1;
                //Calculating the nearest neighbhour for the new sampled point.
                for(unsigned j(0); j<maxSize; j++){
                    if (i==j)
                        continue;
                    float dist = allPoints[i].distance(sampledPoints[j]);
                    if(dist < minDistance){
                        dist = minDistance;
                        minDistIndex = j;
                    }

                }
                nearestNeighbhourDistance[i] = minDistance;
                nearestNeighbhourIndex[i] = minDistIndex;
                sampleSum  += std::log(nearestNeighbhourDistance[i]);
            }
            trueMeanValue = sampleSum/allPoints.size();
        }
        return trueMeanValue;
    }
*/
};

template <class templatePoint>
class Arm2DMutualInformation: public Arm<templatePoint>{
    std::vector<templatePoint> allPoints;
    std::vector<templatePoint> sampledPoints;
    std::vector<unsigned long> indices;
    std::vector<unsigned long> *shuffledRows;
    unsigned long maxSize;

    ArmEntropyContinuous<templatePoint> *Arm11, *Arm10, *Arm01, *Arm1x, *Arm1y;
    long p00, p10, p01, p11 ;


    using Arm<templatePoint>::numberOfPulls;
    using Arm<templatePoint>::sumOfPulls;
    using Arm<templatePoint>::sumOfSquaresOfPulls;
    using Arm<templatePoint>::estimateOfMean;
    using Arm<templatePoint>::estimateOfSecondMoment;
    using Arm<templatePoint>::upperConfidenceBound;
    using Arm<templatePoint>::lowerConfidenceBound;
    using Arm<templatePoint>::localSigma;
    using Arm<templatePoint>::trueMeanValue;

    Arm2DMutualInformation(unsigned long id, std::vector<templatePoint> &allPoints_, std::vector<unsigned long> &indices_){
        indices = indices_;
        allPoints = allPoints_;
        maxSize = allPoints_.size();
        std::vector<unsigned long> shuffledRows_(maxSize);
        std::iota(shuffledRows_.begin(), shuffledRows_.end(), 0);

        std::random_device rd;
        std::mt19937 g(9);
        std::shuffle(shuffledRows_.begin(), shuffledRows_.end(), g);
        shuffledRows = &shuffledRows_;
        ArmEntropyContinuous<templatePoint> Arm11_(5*id, maxSize);
        ArmEntropyContinuous<templatePoint> Arm10_(5*id, maxSize);
        ArmEntropyContinuous<templatePoint> Arm01_(5*id, maxSize);
        ArmEntropyContinuous<templatePoint> Arm1x_(5*id, maxSize);
        ArmEntropyContinuous<templatePoint> Arm1y_(5*id, maxSize);
        Arm11 = &Arm11_;
        Arm10 = &Arm10_;
        Arm01 = &Arm01_;
        Arm1x = &Arm1x_;
        Arm1y = &Arm1y_;

        trueMeanValue = INFINITY;
    }

    templatePoint samplePoint(){
        std::vector<float> sampledVec = allPoints[shuffledRows->at(numberOfPulls)].point;
        std::vector<float> v;
        for(unsigned long i(0); i< indices.size(); i++){
            v.push_back(sampledVec[indices[i]]);
        }
        numberOfPulls++;
        return templatePoint(v);
    }

    using Arm<templatePoint>::pullArm;
    virtual std::pair<float, float> pullArm(float globalSigma, unsigned long globalNumberOfPulls,
                                            float logDeltaInverse, bool update, unsigned sampleSize,
                                            float LCBofSecondBestArm) {
        std::vector<templatePoint> pointsVec11, pointsVec01, pointsVec10, pointsVec1x, pointsVec1y ;
        std::pair<float, float> sample11, sample01, sample10, sample1x, sample1y;

        for(unsigned i(0); i<sampleSize; i++){
            templatePoint p = samplePoint();
            float x = p.point[0];
            float y = p.point[1];
            std::vector<int> tmpx(1,x);
            std::vector<int> tmpy(1,y);
            if( (x==0) and (y==0)){
                p00 ++;
            }
            else if( (y!=0) and (x==0)){
                p01 ++;
                pointsVec01.push_back(templatePoint(tmpy));//H(X,Y)
                pointsVec1y.push_back(templatePoint(tmpy));//H(Y)
            }
            else if( (x!=0) and (y==0)){
                p10 ++;
                pointsVec10.push_back(templatePoint(tmpx));//H(X,Y)
                pointsVec1x.push_back(templatePoint(tmpx));//H(X)
            }
            else{
                p11 ++;
                pointsVec11.push_back(p);//H(X,Y)
                pointsVec1x.push_back(templatePoint(tmpx));//H(X)
                pointsVec1y.push_back(templatePoint(tmpy));//H(Y)
            }
        }

        sample11 = Arm11->pullArm(pointsVec11, globalSigma, logDeltaInverse, update,  LCBofSecondBestArm);
        sample10 = Arm10->pullArm(pointsVec10, globalSigma, logDeltaInverse, update,  LCBofSecondBestArm);
        sample01 = Arm01->pullArm(pointsVec01, globalSigma, logDeltaInverse, update,  LCBofSecondBestArm);
        sample1x = Arm1x->pullArm(pointsVec1x, globalSigma, logDeltaInverse, update,  LCBofSecondBestArm);
        sample1y = Arm1y->pullArm(pointsVec1y, globalSigma, logDeltaInverse, update,  LCBofSecondBestArm);

        if (update){
            updateConfidenceIntervals(globalSigma, logDeltaInverse);
        };
    }
        using Arm<templatePoint>::updateConfidenceIntervals;
        void updateConfidenceIntervals(float globalSigma, float logDeltaInverse){
            float compositeSigma, intervalWidth;
            estimateOfMean = -((p00+0.0)/numberOfPulls)*std::log((p00+0.0)/numberOfPulls)
                             +((p01+0.0)/numberOfPulls)*Arm01->estimateOfMean
                             +((p10+0.0)/numberOfPulls)*Arm10->estimateOfMean
                             +((p11+0.0)/numberOfPulls)*Arm11->estimateOfMean
                             -((p00+p01+0.0)/numberOfPulls)*std::log((p00+p01+0.0)/numberOfPulls)
                             +((p11+p10+0.0)/numberOfPulls)*Arm1x->estimateOfMean
                             -((p00+p10+0.0)/numberOfPulls)*std::log((p00+p10+0.0)/numberOfPulls)
                             +((p11+p01+0.0)/numberOfPulls)*Arm1y->estimateOfMean;

            upperConfidenceBound =
                    -((p00+0.0)/numberOfPulls)*std::log((p00+0.0)/numberOfPulls)*(1+std::pow(numberOfPulls,-0.5))
                    +((p01+0.0)/numberOfPulls)*Arm01->upperConfidenceBound
                    +((p10+0.0)/numberOfPulls)*Arm10->upperConfidenceBound
                    +((p11+0.0)/numberOfPulls)*Arm11->upperConfidenceBound
                    -((p00+p01+0.0)/numberOfPulls)*std::log((p00+p01+0.0)/numberOfPulls)*(1+std::pow(numberOfPulls,-0.5))
                    +((p11+p10+0.0)/numberOfPulls)*Arm1x->upperConfidenceBound
                    -((p00+p10+0.0)/numberOfPulls)*std::log((p00+p10+0.0)/numberOfPulls)*(1+std::pow(numberOfPulls,-0.5))
                    +((p11+p01+0.0)/numberOfPulls)*Arm1y->upperConfidenceBound;

            lowerConfidenceBound =
                    -((p00+0.0)/numberOfPulls)*std::log((p00+0.00)/numberOfPulls)*(1-std::pow(numberOfPulls,-0.5))+
                    +((p01+0.0)/numberOfPulls)*Arm01->lowerConfidenceBound
                    +((p10+0.0)/numberOfPulls)*Arm10->lowerConfidenceBound
                    +((p11+0.0)/numberOfPulls)*Arm11->lowerConfidenceBound
                    -((p00+p01+0.0)/numberOfPulls)*std::log((p00+p01+0.0)/numberOfPulls)*(1-std::pow(numberOfPulls,-0.5))
                    +((p11+p10+0.0)/numberOfPulls)*Arm1x->lowerConfidenceBound
                    -((p00+p10+0.0)/numberOfPulls)*std::log((p00+p10+0.0)/numberOfPulls)*(1-std::pow(numberOfPulls,-0.5))
                    +((p11+p01+0.0)/numberOfPulls)*Arm1y->lowerConfidenceBound;

        }


    };
#endif //COMBINATORIAL_MAB_ARMS_H