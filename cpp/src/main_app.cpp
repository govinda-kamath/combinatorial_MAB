#include"Test.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <cassert>
#include <algorithm>
#include <string>
#include <cstdlib>
#include <queue>
#include <utility>


class Point {
public:
    std::vector <float> point;
    Point(std::vector <float> p){
        point = p;
    }
//    ~Point () {
//        delete &point;
//}

    virtual float distance(const Point& p1){}
    virtual float sampledDistance(const Point& p1){}
};

class SquaredEuclideanPoint : public Point{

public:
    SquaredEuclideanPoint(std::vector<float> p) : Point(p){}

    float distance(const SquaredEuclideanPoint& p1){

        assert(("Sizes do not match", point.size() == p1.point.size()));

        float result(0);

        std::vector<float>::const_iterator pIt = point.begin();
        std::vector<float>::const_iterator p1It = p1.point.begin();
        for (; p1It != p1.point.end() && pIt  != point.end(); ++p1It, ++pIt){
            result += (*p1It-*pIt)*(*p1It-*pIt);
        }
        return result;
    }

    float sampledDistance(const SquaredEuclideanPoint& p1){
        assert(("Sizes do not match", point.size() == p1.point.size()));

        unsigned vecSize = p1.point.size();
        float result;
        unsigned randomCoOrdinate;
        randomCoOrdinate = std::rand() % vecSize;
        return (point[randomCoOrdinate] - p1.point[randomCoOrdinate])
               *(point[randomCoOrdinate] - p1.point[randomCoOrdinate]);
    }
};

template <class templatePoint>
class Arm {
public:
    int numberOfPulls;
    float sumOfPulls;
    float upperConfidenceBound;
    float lowerConfidenceBound;
    float estimateOfMean;
    float estimateOfSecondMoment;
    float SumOfSquaresOfPulls;
    templatePoint * point;

    Arm(){
        numberOfPulls = 0;
        sumOfPulls = 0.0;
        upperConfidenceBound = INFINITY;
        lowerConfidenceBound = -INFINITY;
        SumOfSquaresOfPulls = 0.0;
        estimateOfMean = NAN;
        estimateOfSecondMoment = NAN;
    }

    Arm(templatePoint p){
        Arm();
        point = new templatePoint(p.point);
    }

//    ~Arm(){
//        delete &numberOfPulls;
//        delete &sumOfPulls;
//        delete &upperConfidenceBound;
//        delete &lowerConfidenceBound;
//        delete &estimateOfMean;
//        delete &estimateOfSecondMoment;
//        delete &SumOfSquaresOfPulls;
//        delete &point;
//    }

    void printArm(){
        std::cout << "Number of pulls" << numberOfPulls
                                       << std::endl;
    }

    friend bool operator< (const Arm& l, const Arm& r)
        {
            return l.lowerConfidenceBound < r.lowerConfidenceBound;
        }

    float updateConfidenceIntervals(float globalSigma, float logDeltaInverse){

        float localSigma, intervalWidth;
        localSigma = globalSigma; //Todo: update sigma to new local value
        intervalWidth = std::sqrt((localSigma * localSigma * logDeltaInverse)/numberOfPulls);
        upperConfidenceBound = estimateOfMean + intervalWidth;
        lowerConfidenceBound = std::min((float)0.0, estimateOfMean - intervalWidth);
    }

    float pullArm(const templatePoint &p1, float globalSigma,
                            float logDeltaInverse, bool update = true) {
        float sample;

        sample = point->sampledDistance(p1);

        numberOfPulls ++;
        sumOfPulls += sample;
        estimateOfMean = sumOfPulls/numberOfPulls;
        SumOfSquaresOfPulls += sample*sample;
        estimateOfSecondMoment = SumOfSquaresOfPulls/numberOfPulls;

        if (update)
            updateConfidenceIntervals(globalSigma,logDeltaInverse);

        return sample;
    }


    float trueMean(const templatePoint &p1){
        return point->distance(p1);
    }


};


template <class templatePoint>
class ArmKNN : public Arm<templatePoint>{
public:
    templatePoint *fixedPoint;

    ArmKNN(templatePoint p) : Arm<templatePoint>(p) {}

    ArmKNN(templatePoint p, templatePoint fixPoint) : Arm<templatePoint>(p) {
        fixedPoint = new templatePoint(fixPoint.point);
    }

    using Arm<templatePoint>::pullArm;
    float pullArm(float globalSigma, float logDeltaInverse, bool update = true){
        pullArm(*fixedPoint, globalSigma, logDeltaInverse, update);
    }

};

template <class templatePoint>
bool armsLessThan(const Arm<templatePoint>& l, const Arm<templatePoint>& r)
{
    return l.lowerConfidenceBound < r.lowerConfidenceBound;
}

template <class templateArm>
class UCB{
public:
    unsigned numberOfArms;
    std::vector<templateArm> armsContainer;
    std::priority_queue<templateArm, std::vector<templateArm> > arms;
    float logDeltaInverse;

    float globalSigma;
    float globalNumberOfPulls;
    float globalSumOfPulls;
    float globalSumOfSquaresOfPulls;


    UCB(std::vector<templateArm> armsVec, float delta){

        armsContainer = armsVec;
        numberOfArms = armsContainer.size();

        logDeltaInverse = std::log(1/delta);

        globalSigma = NAN;
        globalNumberOfPulls = 0;
        globalSumOfPulls = 0;
        globalSumOfSquaresOfPulls = 0;

    }

    void initialise(int numberOfInitialPulls = 100){

        for (typename std::vector<templateArm>::const_iterator it = armsContainer.begin();
             it != armsContainer.end(); ++it){
            for (unsigned i = 0; i < numberOfInitialPulls; i++) {
                float observedSample;
                observedSample = it->pullArm(NAN, NAN, false);
                globalSumOfPulls += observedSample;
                globalSumOfSquaresOfPulls += observedSample * observedSample;
            }
        }

        globalNumberOfPulls = numberOfInitialPulls*armsContainer.size();
        globalSigma = (globalSumOfSquaresOfPulls - std::pow(globalSumOfPulls,2))/globalNumberOfPulls;

        for (typename std::vector<templateArm>::const_iterator it = armsContainer.begin();
             it != armsContainer.end(); ++it){
            it->updateConfidenceIntervals(globalSigma, logDeltaInverse);
            arms.push(*it);
        }
    }


};



int main(int argc, char *argv[]){

    float testArray[4] = {0, 1, 3, 5};
    std::vector<float> testVector ;
    for (int index=0; index<4; index++) {
        testVector.push_back(index);
        std::cout << index <<std::endl;
    }
    SquaredEuclideanPoint testPoint(testVector);
    return 0;
}
