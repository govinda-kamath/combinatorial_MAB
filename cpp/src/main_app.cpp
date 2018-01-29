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
    ~Point () {
        delete &point;
}

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


class Arm {
public:
    int numberOfPulls;
    float sumOfPulls;
    float upperConfidenceBound;
    float lowerConfidenceBound;
    float estimateOfMean;
    float estimateOfSecondMoment;
    float SumOfSquaresOfPulls;
    Point * point;

    Arm(){
        numberOfPulls = 0;
        sumOfPulls = 0.0;
        upperConfidenceBound = INFINITY;
        lowerConfidenceBound = -INFINITY;
        SumOfSquaresOfPulls = 0.0;
        estimateOfMean = NAN;
        estimateOfSecondMoment = NAN;
    }

    Arm(SquaredEuclideanPoint p){
        Arm();
        point = new SquaredEuclideanPoint(p.point);
    }

    ~Arm(){
        delete &numberOfPulls;
        delete &sumOfPulls;
        delete &upperConfidenceBound;
        delete &lowerConfidenceBound;
        delete &estimateOfMean;
        delete &estimateOfSecondMoment;
        delete &SumOfSquaresOfPulls;
        delete &point;
    }

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

    float pullArm(const SquaredEuclideanPoint &p1, float globalSigma,
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


    float trueMean(const SquaredEuclideanPoint &p1){
        return point->distance(p1);
    }


};

template <class T>

class ArmKNN : public Arm{
public:
    T *fixedPoint;

    ArmKNN(T p) : Arm(p) {}

    ArmKNN(T p, T fixPoint) : Arm(p) {
        fixedPoint = new T(fixPoint.point);
    }

    using Arm::pullArm;
    float pullArm(float globalSigma, float logDeltaInverse, bool update = true){
        pullArm(*fixedPoint, globalSigma, logDeltaInverse, update);
    }

};

bool armsLessThan(const Arm& l, const Arm& r)
{
    return l.lowerConfidenceBound < r.lowerConfidenceBound;
}

class UCB{
public:
    unsigned numberOfArms;
    std::vector<Arm> armsContainer;
    std::priority_queue<Arm, std::vector<Arm> > arms;
    float logDeltaInverse;

    float globalSigma;
    float globalNumberOfPulls;
    float globalSumOfPulls;
    float globalSumOfSquaresOfPulls;


    UCB(std::vector<Arm> armsContainer, float delta){

        numberOfArms = armsContainer.size();

//        for (std::vector<Arm>::const_iterator it = armsContainer.begin(); it != armsContainer.end(); ++it){
////            it->pullArm()
//        }

        for (std::vector<Arm>::const_iterator it = armsContainer.begin(); it != armsContainer.end(); ++it)
            arms.push(*it);

        logDeltaInverse = std::log(1/delta);

        globalSigma = NAN;
        globalNumberOfPulls = 0;
        globalSumOfPulls = 0;
        globalSumOfSquaresOfPulls = 0;
    }



};



int main(int argc, char *argv[]){
   Test s("Joe");
   s.display();
   return 0;
}
