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
#include <math.h>       /* tgamma */
#include <boost/math/special_functions/digamma.hpp>

//#include <nanoflann.hpp>
#define chi_d 0.0 //the constant multiplying the variance

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
class ArmEntropyContinuousHelper: public Arm<templatePoint>{
public:
    std::vector<templatePoint> sampledPoints;
    std::vector<float> nearestNeighbhourDistance;
    std::vector<unsigned long> nearestNeighbhourIndex;
    unsigned long maxSize;
    float constant, intervalWidth;
    unsigned dim;

    using Arm<templatePoint>::numberOfPulls;
    using Arm<templatePoint>::sumOfPulls; //summation ln(pi)
    using Arm<templatePoint>::sumOfSquaresOfPulls;//summation ln^2(pi)
    using Arm<templatePoint>::estimateOfMean;
    using Arm<templatePoint>::estimateOfSecondMoment;
    using Arm<templatePoint>::upperConfidenceBound;
    using Arm<templatePoint>::lowerConfidenceBound;
    using Arm<templatePoint>::localSigma;
    using Arm<templatePoint>::trueMeanValue;
    using Arm<templatePoint>::id;


    ArmEntropyContinuousHelper() : Arm<templatePoint>(){
        estimateOfMean = 0;
        upperConfidenceBound = 0;
        lowerConfidenceBound = 0;
        sumOfPulls = 0;
        sumOfSquaresOfPulls = 0;
    }

    ArmEntropyContinuousHelper(unsigned long armNumber) : Arm<templatePoint>(){
        id = armNumber;
    }

    void initialize(long id_, unsigned long maxSize_, float dimension_){
        id = id_;
        maxSize = maxSize_;
        for(unsigned long i(0); i< maxSize; i++) {
            nearestNeighbhourIndex.push_back(-1);
            nearestNeighbhourDistance.push_back(INFINITY);
        }
        dim = (unsigned) dimension_;
        constant = (float)std::log( std::pow(M_PI,dimension_/2)/std::tgamma(dimension_/2+1) );

    }

    using Arm<templatePoint>::pullArm;
    virtual std::pair<float, float> pullArm(std::vector<templatePoint> samplePoints, float globalSigma,
                                            float logDeltaInverse, bool update, float LCBofSecondBestArm) {
        float sampleSum(0), sampleSquareSum(0);
        for (unsigned i(0); i < samplePoints.size(); i++) {

            templatePoint p = samplePoints[i];
//            if(dim==1)
//                std::cout << p.point[0] << std::endl;
            float minDistance(INFINITY);
            unsigned long minDistIndex(-1);
            //Calculating the nearest neighbhour for the new sampled point.
            for (unsigned j(0); j < numberOfPulls; j++) {
                float dist = std::sqrt(p.distance(sampledPoints[j]));
                if(dist==0)
                    continue;
                if (dist < minDistance) {
                    minDistance = dist;
                    minDistIndex = j;
                }

                if (dist < nearestNeighbhourDistance[j]) {
                    if(nearestNeighbhourDistance[j] != INFINITY){
                        sampleSum -= std::log(nearestNeighbhourDistance[j]);
                        sampleSquareSum -= std::pow(std::log(nearestNeighbhourDistance[j]), 2);
                    }
                    sampleSum += std::log(dist);
                    sampleSquareSum += std::pow(std::log(dist), 2);

                    nearestNeighbhourDistance[j] = dist;
                    nearestNeighbhourIndex[j] = numberOfPulls;
                }
            }
            if (minDistance != INFINITY) {
                sampleSum += std::log(minDistance);
                sampleSquareSum += std::pow(std::log(minDistance), 2);
                nearestNeighbhourDistance[numberOfPulls] = minDistance;
                nearestNeighbhourIndex[numberOfPulls] = minDistIndex;
            }
            sampledPoints.push_back(p);
            numberOfPulls++;
        }

        if (sampleSum != 0) {
            sumOfPulls += sampleSum;
            sumOfSquaresOfPulls += sampleSquareSum;
        }
        estimateOfMean = dim*sumOfPulls/numberOfPulls;
        estimateOfSecondMoment = dim*dim*sumOfSquaresOfPulls/numberOfPulls;
        float tmpSum = 0;
//            for( int jj(0); jj < numberOfPulls; jj++){
//                tmpSum += nearestNeighbhourDistance[jj];
//            }
//        std::cout << "Id " << id << " avg N distance  " << tmpSum/numberOfPulls << " " << numberOfPulls << std::endl;


        if (numberOfPulls>0){
            updateConfidenceIntervals(globalSigma, logDeltaInverse);
        }
        return std::make_pair(sampleSum, sampleSquareSum);
    }


    using Arm<templatePoint>::updateConfidenceIntervals;
    void updateConfidenceIntervals(float globalSigma, float logDeltaInverse){
        float compositeSigma;
        float tmp = std::abs(estimateOfSecondMoment - std::pow(estimateOfMean,2));
        if( numberOfPulls==2)
            compositeSigma = 0;
        else {
            localSigma = std::sqrt(estimateOfSecondMoment - std::pow(estimateOfMean, 2));
            compositeSigma = std::sqrt(localSigma * localSigma + chi_d);
        }
        if (localSigma<0){
            std::cout << "Abort mission!! Fundamental error" <<std::endl;
        }
        float frac = (float) numberOfPulls/maxSize;
        if (frac>=1){
            intervalWidth = 0;
        }
        else{
//            compositeSigma = std::sqrt(localSigma*localSigma*frac +  globalSigma*globalSigma*(1- frac)+ chi_d);
            intervalWidth = std::sqrt((compositeSigma * compositeSigma * logDeltaInverse)/(float)numberOfPulls);
        }
        estimateOfMean += boost::math::digamma(numberOfPulls)-boost::math::digamma(1)+constant;
        upperConfidenceBound = estimateOfMean + intervalWidth;
        lowerConfidenceBound = estimateOfMean - intervalWidth;
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
class Arm1DEntropy: public Arm<templatePoint>{
public:
    std::vector<templatePoint> allPoints;
    unsigned long index;
    std::vector<unsigned long> shuffledRows;
    unsigned long maxSize;

    ArmEntropyContinuousHelper<templatePoint> Arm1;
    using Arm<templatePoint>::numberOfPulls;
    using Arm<templatePoint>::sumOfPulls;
    using Arm<templatePoint>::sumOfSquaresOfPulls;
    using Arm<templatePoint>::estimateOfMean;
    using Arm<templatePoint>::estimateOfSecondMoment;
    using Arm<templatePoint>::upperConfidenceBound;
    using Arm<templatePoint>::lowerConfidenceBound;
    using Arm<templatePoint>::localSigma;
    using Arm<templatePoint>::trueMeanValue;
    using Arm<templatePoint>::id;

    Arm1DEntropy(long id_) {
        id = id_;
    }

    Arm1DEntropy(long id_, std::vector<templatePoint> &allPoints_,
                           unsigned long index_, std::vector<unsigned long> shuffledRows_) {
        id = id_;
        index = index_;
        allPoints = allPoints_;
        maxSize = allPoints_.size();
        shuffledRows = shuffledRows_;
        Arm1.initialize(id, maxSize, 1);
    }

    templatePoint samplePoint(){
        std::vector<float> v = {allPoints[shuffledRows[numberOfPulls]].point[index]};
        numberOfPulls++;
        return templatePoint(v);
    }

    using Arm<templatePoint>::pullArm;
    virtual std::pair<float, float> pullArm(float globalSigma, unsigned long globalNumberOfPulls,
                                            float logDeltaInverse, bool update, unsigned sampleSize,
                                            float LCBofSecondBestArm) {
        std::vector<templatePoint> pointsVec;
        std::pair<float, float> sample;
        for(unsigned i(0); i<sampleSize; i++){
            templatePoint p = samplePoint();
            pointsVec.push_back(p);
        }

        sample = Arm1.pullArm(pointsVec, globalSigma, logDeltaInverse, update,  LCBofSecondBestArm);
        estimateOfMean = Arm1.estimateOfMean;
        lowerConfidenceBound = Arm1.lowerConfidenceBound;
        upperConfidenceBound = Arm1.upperConfidenceBound;

        std::cout << id << " " << numberOfPulls
                << "\t" << Arm1.estimateOfMean
                << " " << Arm1.lowerConfidenceBound
                << " " << Arm1.upperConfidenceBound
                << std::endl;
        return sample;
    }

};

template <class templatePoint>
class Arm2DMutualInformation: public Arm<SquaredEuclideanPoint>{
public:
    std::vector<templatePoint> allPoints;
    std::vector<unsigned long> indices;
    unsigned long maxSize;
    float intervalWidth;

//    ArmEntropyContinuousHelper<templatePoint> *Arm11, *Arm10, *Arm01, *Arm1x, *Arm1y;
    ArmEntropyContinuousHelper<SquaredEuclideanPoint> Arm11;
    ArmEntropyContinuousHelper<SquaredEuclideanPoint> Arm10;
    ArmEntropyContinuousHelper<SquaredEuclideanPoint> Arm01;
    ArmEntropyContinuousHelper<SquaredEuclideanPoint> Arm1x;
    ArmEntropyContinuousHelper<SquaredEuclideanPoint> Arm1y;
    long p00, p10, p01, p11;


    using Arm<SquaredEuclideanPoint>::numberOfPulls;
    using Arm<SquaredEuclideanPoint>::sumOfPulls;
    using Arm<SquaredEuclideanPoint>::sumOfSquaresOfPulls;
    using Arm<SquaredEuclideanPoint>::estimateOfMean;
    using Arm<SquaredEuclideanPoint>::estimateOfSecondMoment;
    using Arm<SquaredEuclideanPoint>::upperConfidenceBound;
    using Arm<SquaredEuclideanPoint>::lowerConfidenceBound;
    using Arm<SquaredEuclideanPoint>::localSigma;
    using Arm<SquaredEuclideanPoint>::trueMeanValue;
    using Arm<SquaredEuclideanPoint>::id;

    Arm2DMutualInformation(unsigned long id_) {
        id = id_;
    }

    Arm2DMutualInformation(unsigned long id_, std::vector<templatePoint> &allPoints_,
                           std::vector<unsigned long> indices_) {
        id = id_;
        indices = indices_;
        allPoints = allPoints_;
        maxSize = allPoints_.size();
        Arm11.initialize(10*id+0, maxSize, 2);
        Arm10.initialize(10*id+1, maxSize, 1);
        Arm01.initialize(10*id+2, maxSize, 1);
        Arm1x.initialize(10*id+3, maxSize, 1);
        Arm1y.initialize(10*id+4, maxSize, 1);
        p00 = 0;
        p10 = 0;
        p01 = 0;
        p11 = 0;
    }

    SquaredEuclideanPoint samplePoint(){
//        std::unordered_map <unsigned long, float> sampledVec = allPoints[shuffledRows[numberOfPulls]].sparsePoint;
        std::vector<float> sampledVec = allPoints[numberOfPulls].point;
        std::vector<float> v;
        for(unsigned long i(0); i< indices.size(); i++){
            v.push_back(sampledVec[indices[i]]);
        }
        numberOfPulls++;
        return SquaredEuclideanPoint(v);
    }

    float xlogx(float x){
        if (x==0)
            return 0;
        else
            return x*std::log(x);
    }

    using Arm<SquaredEuclideanPoint>::pullArm;
    virtual std::pair<float, float> pullArm(float globalSigma, unsigned long globalNumberOfPulls,
                                            float logDeltaInverse, bool update, unsigned sampleSize,
                                            float LCBofSecondBestArm) {
        std::vector<SquaredEuclideanPoint> pointsVec11, pointsVec01, pointsVec10, pointsVec1x, pointsVec1y ;
        std::pair<float, float> sample11, sample01, sample10, sample1x, sample1y;
        unsigned Lp00 = 0;
        unsigned Lp10 = 0;
        unsigned Lp01 = 0;
        unsigned Lp11 = 0;
        for(unsigned i(0); i<sampleSize; i++){
            SquaredEuclideanPoint p = samplePoint();
            float x = p.point[0];
            float y = p.point[1];
//            std::cout << x << " " << y << std::endl;
            std::vector<float> tmpx = {x};
            std::vector<float> tmpy = {y};
            if( (x==0) and (y==0)){
                p00 ++;
                Lp00 ++;
            }
            else if((y!=0) and (x==0)){
                p01 ++;
                Lp01 ++;
                pointsVec01.push_back(SquaredEuclideanPoint(tmpy));//H(X,Y)
                pointsVec1y.push_back(SquaredEuclideanPoint(tmpy));//H(Y)
            }
            else if((x!=0) and (y==0)){
                p10 ++;
                Lp10 ++;
                pointsVec10.push_back(SquaredEuclideanPoint(tmpx));//H(X,Y)
                pointsVec1x.push_back(SquaredEuclideanPoint(tmpx));//H(X)
            }
            else{
                p11 ++;
                Lp11 ++;
                pointsVec11.push_back(p);//H(X,Y)
                pointsVec1x.push_back(SquaredEuclideanPoint(tmpx));//H(X)
                pointsVec1y.push_back(SquaredEuclideanPoint(tmpy));//H(Y)
            }
        }

        sample11 = Arm11.pullArm(pointsVec11, globalSigma, logDeltaInverse, update,  LCBofSecondBestArm);
        if(p10>=2)
            sample10 = Arm10.pullArm(pointsVec10, globalSigma, logDeltaInverse, update,  LCBofSecondBestArm);
        if(p01>=2)
            sample01 = Arm01.pullArm(pointsVec01, globalSigma, logDeltaInverse, update,  LCBofSecondBestArm);
        sample1x = Arm1x.pullArm(pointsVec1x, globalSigma, logDeltaInverse, update,  LCBofSecondBestArm);
        sample1y = Arm1y.pullArm(pointsVec1y, globalSigma, logDeltaInverse, update,  LCBofSecondBestArm);


        if(p01==1) {
            Arm01.lowerConfidenceBound = Arm01.estimateOfMean;
            Arm01.upperConfidenceBound = Arm01.estimateOfMean;
        }
        if(p10==1) {
            Arm10.lowerConfidenceBound = Arm01.estimateOfMean;
            Arm10.upperConfidenceBound = Arm01.estimateOfMean;
        }

        if(Arm10.lowerConfidenceBound==NAN){
            std::cout << " NAN here! " << std::endl;
        }

        if(Arm01.lowerConfidenceBound==NAN){
            std::cout << " NAN here! " << std::endl;
        }

        float sumOfPullsPrev = sumOfPulls;
        float sumOfSquaresOfPullsPrev = sumOfSquaresOfPulls;
        updateConfidenceIntervals(globalSigma, 0, logDeltaInverse);
        float first = sumOfPulls - sumOfPullsPrev;
        float second = sumOfSquaresOfPulls - sumOfSquaresOfPullsPrev;

        std::cout << "ID " << id
                    << "\n11\t" <<   Arm11.lowerConfidenceBound
                    << " " <<      Arm11.estimateOfMean
                    << " " <<      Arm11.upperConfidenceBound
                    << "\n01\t" <<   Arm01.lowerConfidenceBound
                    << " " <<      Arm01.estimateOfMean
                    << " " <<      Arm01.upperConfidenceBound
                    << "\n10\t" <<   Arm10.lowerConfidenceBound
                    << " " <<      Arm10.estimateOfMean
                    << " " <<      Arm10.upperConfidenceBound
                    << "\nx\t" <<   Arm1x.lowerConfidenceBound
                    << " " <<      Arm1x.estimateOfMean
                    << " " <<      Arm1x.upperConfidenceBound
                    << "\ny\t" <<   Arm1y.lowerConfidenceBound
                    << " " <<      Arm1y.estimateOfMean
                    << " " <<      Arm1y.upperConfidenceBound
                   << "\nF\t" <<   lowerConfidenceBound
                   << " " <<      estimateOfMean
                   << " " <<      upperConfidenceBound
                           <<std::endl;

        return std::make_pair(first, second);
    }
    using Arm<SquaredEuclideanPoint>::updateConfidenceIntervals;
    void updateConfidenceIntervals(float globalSigma, unsigned long globalNumberOfPulls, float logDeltaInverse) {
        float compositeSigma;
        estimateOfMean = -(-xlogx((p00 + 0.0) / numberOfPulls)
                           + ((p01 + 0.0) / numberOfPulls) * Arm01.estimateOfMean
                           + ((p10 + 0.0) / numberOfPulls) * Arm10.estimateOfMean
                           + ((p11 + 0.0) / numberOfPulls) * Arm11.estimateOfMean
        ) // H(X,Y) till here
                         - (xlogx((p00 + p01 + 0.0) / numberOfPulls))
                         + ((p11 + p10 + 0.0) / numberOfPulls) * Arm1x.estimateOfMean
                         - (xlogx((p00 + p10 + 0.0) / numberOfPulls))
                         + ((p11 + p01 + 0.0) / numberOfPulls) * Arm1y.estimateOfMean;

        upperConfidenceBound =
                -(-xlogx((p00 + 0.0) / numberOfPulls) *
                  (1 - std::pow(numberOfPulls, -0.5) * std::sqrt((p00 + 0.0) / numberOfPulls))
                  + ((p01 + 0.0) / numberOfPulls) * Arm01.lowerConfidenceBound
                  + ((p10 + 0.0) / numberOfPulls) * Arm10.lowerConfidenceBound
                  + ((p11 + 0.0) / numberOfPulls) * Arm11.lowerConfidenceBound)
                - (xlogx((p00 + p01 + 0.0) / numberOfPulls)) *
                  (1 + std::pow(numberOfPulls, -0.5) * std::sqrt((p00 + p01 + 0.0) / numberOfPulls))
                + ((p11 + p10 + 0.0) / numberOfPulls) * Arm1x.upperConfidenceBound
                - (xlogx((p00 + p10 + 0.0) / numberOfPulls)) *
                  (1 + std::pow(numberOfPulls, -0.5) * std::sqrt((p00 + p10 + 0.0) / numberOfPulls))
                + ((p11 + p01 + 0.0) / numberOfPulls) * Arm1y.upperConfidenceBound;

        lowerConfidenceBound =
                -(-xlogx((p00 + 0.0) / numberOfPulls) * (1 + std::pow(numberOfPulls, -0.5))
                  + ((p01 + 0.0) / numberOfPulls) * Arm01.upperConfidenceBound
                  + ((p10 + 0.0) / numberOfPulls) * Arm10.upperConfidenceBound
                  + ((p11 + 0.0) / numberOfPulls) * Arm11.upperConfidenceBound)
                - (xlogx((p00 + p01 + 0.0) / numberOfPulls)) *
                  (1 - std::pow(numberOfPulls, -0.5) * std::sqrt((p00 + p01 + 0.0) / numberOfPulls))
                + ((p11 + p10 + 0.0) / numberOfPulls) * Arm1x.lowerConfidenceBound
                - (xlogx((p00 + p10 + 0.0) / numberOfPulls)) *
                  (1 - std::pow(numberOfPulls, -0.5) * std::sqrt((p00 + p10 + 0.0) / numberOfPulls))
                + ((p11 + p01 + 0.0) / numberOfPulls) * Arm1y.lowerConfidenceBound;

        intervalWidth = std::max(estimateOfMean - lowerConfidenceBound, upperConfidenceBound - estimateOfMean);
        estimateOfSecondMoment = estimateOfMean + intervalWidth;
        sumOfPulls = estimateOfMean * numberOfPulls;
        sumOfSquaresOfPulls = estimateOfSecondMoment * numberOfPulls;

        // Because UCB is picking the minimum.
        estimateOfMean *= -1;
        float tmp = upperConfidenceBound;
        upperConfidenceBound = -1 * lowerConfidenceBound;
        lowerConfidenceBound = -1 * tmp;

    }
};
#endif //COMBINATORIAL_MAB_ARMS_H