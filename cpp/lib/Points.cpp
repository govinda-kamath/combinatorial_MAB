//
// Created by Vivek Kumar Bagaria on 1/31/18.
//
//

#include "Points.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <thread>
#include <cassert>

//BasePoint
BasePoint::BasePoint(){
    dimensionPointer = 0;
}
float BasePoint::distance(const BasePoint &p1) const { return  -1; }
std::pair<float, float> BasePoint::sampleDistance(const BasePoint &p1) const { return  std::make_pair(-1.0, -1.0); }
unsigned long BasePoint::getVecSize() const {
    return vecSize;
}

/*Point*/
Point::Point(): BasePoint(){}
Point::Point(const std::vector <float> &p): Point() {
    point = p;
    vecSize = p.size();
}
//Sparse point
SparsePoint::SparsePoint(): BasePoint(){}
SparsePoint::SparsePoint(const std::unordered_map <unsigned long, float>  &sp, unsigned long d): SparsePoint(){
    sparsePoint = sp;
    vecSize = d;
    for(auto kv : sp) {
        keys.push_back(kv.first);
    }
}

/*SquaredEuclideanPoint*/
SquaredEuclideanPoint::SquaredEuclideanPoint(const std::vector<float> &p): Point(p){}
/*Computes the exact distance between two points.
 * Used only for debug purposes*/
float SquaredEuclideanPoint::distance(const SquaredEuclideanPoint &p1) const {
    assert(("Sizes do not match", point.size() == p1.point.size()));

    float result(0);

    std::vector<float>::const_iterator pIt = point.begin();
    std::vector<float>::const_iterator p1It = p1.point.begin();
    for (; p1It != p1.point.end() && pIt  != point.end(); ++p1It, ++pIt){
        result += (*p1It-*pIt)*(*p1It-*pIt);
    }
    return result;
}
/*Picks a dimension of points randomly and samples the distance
 * that dimension*/
std::pair<float, float> SquaredEuclideanPoint::sampleDistance(const SquaredEuclideanPoint &p1, const unsigned sampleSize)  const {
    assert(("Sizes do notmatch", point.size() == p1.point.size()));
    std::srand(std::time(nullptr));
    float sampleSum = 0;
    float sampleSquareSum = 0;
//    unsigned int dimensionPointer = std::rand(); // Todo: Bad code
//    std::cout << dimensionPointer << std::endl;
    for(unsigned i(0); i< sampleSize; i++) {
        unsigned long coordinate = std::rand() % getVecSize();
        float value = (point[coordinate] - p1.point[coordinate])
                     *(point[coordinate] - p1.point[coordinate]);
        sampleSum += value;
        sampleSquareSum += value*value;
//        dimensionPointer ++;
    }
    return std::make_pair(sampleSum, sampleSquareSum);
}

/* L1Point */
L1Point::L1Point(const std::vector<float> &p): Point(p){}
/*Computes the exact distance between two points.
 * Used only for debug purposes*/
float L1Point::distance(const L1Point &p1) const {
    assert(("Sizes do not match", point.size() == p1.point.size()));
    float result(0);

    std::vector<float>::const_iterator pIt = point.begin();
    std::vector<float>::const_iterator p1It = p1.point.begin();
    unsigned i(0);
    for (; p1It != p1.point.end() && pIt  != point.end(); ++p1It, ++pIt){
        i++;
        result += std::abs(*p1It-*pIt);
    }
    return result;
}
/*Picks a dimension of points randomly and samples the distance
 * that L1Point*/
std::pair<float, float> L1Point::sampleDistance(const L1Point &p1, const unsigned sampleSize) {
    assert(("Sizes do not match", point.size() == p1.point.size()));
    float sampleSum = 0;
    float sampleSquareSum = 0;
    for(unsigned i(0); i< sampleSize; i++) {
        unsigned long coordinate = std::rand() % getVecSize();
        float value = std::abs(point[coordinate] - p1.point[coordinate]);
        sampleSum += value;
        sampleSquareSum += value*value;
//        dimensionPointer ++;
    }
    return std::make_pair(sampleSum, sampleSquareSum);
}

/* Sparse L1Point */
SparseL1Point::SparseL1Point(const std::unordered_map <unsigned long, float> &sp, unsigned long d): SparsePoint(sp, d){}

/*Computes the exact distance between two points.
 * Used only for debug purposes*/
float SparseL1Point::distance(const SparseL1Point &p1) const {
//    assert(("Sizes do not match", sparsePoint.size() == p1.sparsePoint.size()));

    float result(0);

    for(unsigned long i(0); i < vecSize; i++){
        // Finding two points
        auto search1 = sparsePoint.find(i);
        auto search2 = p1.sparsePoint.find(i);
        float value(0);
        // If both the points have the index i
        if ( (search1 != sparsePoint.end()) && (search2 != p1.sparsePoint.end()) ){
            value += std::abs(search1->second - search2->second);
        }
        // If only first point has the index i
        else if ( (search1 != sparsePoint.end()) && (search2 == p1.sparsePoint.end() )){
            value += std::abs(search1->second);
        }
        // If only second point has the index i
        else if ( (search1 == sparsePoint.end()) && (search2 != p1.sparsePoint.end() )){
            value += std::abs(search2->second);
        }
        result += value;
    }

    return result;
}

/*Picks a dimension of points randomly and samples the distance
 * that L1Point*/
std::pair<float, float> SparseL1Point::sampleDistance(const SparseL1Point &p1, const unsigned sampleSize ) const {

    float sampleSum = 0;
    float sampleSquareSum = 0;
    float result(0);

    for(unsigned i(0); i< sampleSize; i++) {

        {
            // pick a random index from the current point
            unsigned long randomCoOrdinate = std::rand() % keys.size();
            unsigned long index = keys[randomCoOrdinate];

            // Check if the index exists in the other point
            auto search1 = sparsePoint.find(index);
            auto search2 = p1.sparsePoint.find(index);
            if (search2 != p1.sparsePoint.end())
                result = std::abs(search1->second - search2->second) * keys.size() / getVecSize();
            else
                result = 2 * std::abs(search1->second) * keys.size() / getVecSize();
        }
        sampleSum += result;
        sampleSquareSum += result * result;
        {
            // pick a random index from the current point
            unsigned long randomCoOrdinate = std::rand() % p1.keys.size();
            unsigned long index = p1.keys[randomCoOrdinate];

            // Check if the index exists in the first point
            auto search1 = sparsePoint.find(index);
            auto search2 = p1.sparsePoint.find(index);
            if (search1 != sparsePoint.end())
                result = std::abs(search1->second - search2->second) * p1.keys.size() / getVecSize();
            else
                result = 2 * std::abs(search2->second) * p1.keys.size() / getVecSize();
        }
        sampleSum += result;
        sampleSquareSum += result * result;
    }
    return std::make_pair(sampleSum/2, sampleSquareSum/2);
}