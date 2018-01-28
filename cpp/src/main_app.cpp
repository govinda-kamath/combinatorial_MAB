#include"Test.h"
#include <cmath>
#include <iostream>
#include <vector>

class Point {
    public:
        std::vector <float> point;
    Point(std::vector <float> p){
        point = p;
    }
    ~Point () {
        delete point;
    }

    friend float distance(const Point& p1, const Point& p2){
        return 0;
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
        Point point;

    Arm(Point p){
        numberOfPulls = 0;
        sumOfPulls = 0.0;
        upperConfidenceBound = INFINITY;
        lowerConfidenceBound = -INFINITY;
        estimateOfMean = NAN;
        estimateOfSecondMoment = NAN;
        point = p;
    }

    ~Arm(){
        delete numberOfPulls;
        delete sumOfPulls;
        delete upperConfidenceBound;
        delete lowerConfidenceBound;
        delete estimateOfMean;
        delete estimateOfSecondMoment;
        delete point;
    }

    void printArm(){
        std::cout << "Number of pulls" << numberOfPulls
                                       << std::endl;
    }

    friend bool operator<(const Arm& l, const Arm& r)
    {
        return l.lowerConfidenceBound < r.lowerConfidenceBound;
    }

};

class UCB{
    public:


};



int main(int argc, char *argv[]){
   Test s("Joe");
   s.display();
   return 0;
}
