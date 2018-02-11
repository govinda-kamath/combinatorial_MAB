//
// Created by Vivek Kumar Bagaria on 1/31/18.
//

#include "Points.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <thread>
#include <cassert>

//BasePoint
float BasePoint::distance(const BasePoint& p1) const { return  -1; }
float BasePoint::sampledDistance(const BasePoint& p1) const { return  -1; }
unsigned long BasePoint::getVecSize() const {
    return vecSize;
}

//Point
Point::Point(std::vector <float> p){
    point = p;
    vecSize = p.size();
}
//Sparse point
SparsePoint::SparsePoint(std::unordered_map <unsigned long, int>  sp){
    sparsePoint = sp;
    vecSize = sp.size();
}

// SquaredEuclideanPoint
SquaredEuclideanPoint::SquaredEuclideanPoint(std::vector<float> p) : Point(p){}
/*Computes the exact distance between two points.
 * Used only for debug purposes*/
float SquaredEuclideanPoint::distance(const SquaredEuclideanPoint& p1) const {
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
float SquaredEuclideanPoint::sampledDistance(const SquaredEuclideanPoint& p1) const {
    assert(("Sizes do not match", point.size() == p1.point.size()));
    unsigned vecSize = getVecSize();
    unsigned randomCoOrdinate;
    randomCoOrdinate = std::rand() % vecSize;
    return (point[randomCoOrdinate] - p1.point[randomCoOrdinate])
           *(point[randomCoOrdinate] - p1.point[randomCoOrdinate]);
}

// L1Point
L1Point::L1Point(std::vector<float> p) : Point(p){}
/*Computes the exact distance between two points.
 * Used only for debug purposes*/
float L1Point::distance(const L1Point& p1) const {
    assert(("Sizes do not match", point.size() == p1.point.size()));

    float result(0);

    std::vector<float>::const_iterator pIt = point.begin();
    std::vector<float>::const_iterator p1It = p1.point.begin();
    for (; p1It != p1.point.end() && pIt  != point.end(); ++p1It, ++pIt){
        result += std::abs(*p1It-*pIt);
    }
    return result;
}
/*Picks a dimension of points randomly and samples the distance
 * that L1Point*/
float L1Point::sampledDistance(const L1Point& p1) const {
    assert(("Sizes do not match", point.size() == p1.point.size()));
    unsigned randomCoOrdinate;
    randomCoOrdinate = std::rand() % getVecSize();
    return std::abs(point[randomCoOrdinate] - p1.point[randomCoOrdinate]);

}


SparseL1Point::SparseL1Point(std::unordered_map <unsigned long, int> sp): SparsePoint(sp){
    sparsePoint = sp;
    for(auto kv : sp) {
        keys.push_back(kv.first);
    }
}
/*Computes the exact distance between two points.
 * Used only for debug purposes*/
float SparseL1Point::distance(const SparseL1Point& p1) const {
    assert(("Sizes do not match", sparsePoint.size() == p1.sparsePoint.size()));

    float result(0);

    for(unsigned long i(0); i < sparsePoint.size(); i++){
        auto search1 = sparsePoint.find(i);
        auto search2 = p1.sparsePoint.find(i);

        // If both the points have the index i
        if ( (search1 != sparsePoint.end()) && (search2 != p1.sparsePoint.end() )){
            result += std::abs(search1->second - search2->second);
        }
        // If only first point has the index i
        else if ( (search1 != sparsePoint.end()) && (search2 == p1.sparsePoint.end() )){
            result += std::abs(search1->second);
        }
        // If only second point has the index i
        else if ( (search1 == sparsePoint.end()) && (search2 != p1.sparsePoint.end() )){
            result += std::abs(search2->second);
        }
        //If none of the points have the index, we do not update result.

    }

    return result;
}

/*Picks a dimension of points randomly and samples the distance
 * that L1Point*/
float SparseL1Point::sampledDistance(const SparseL1Point& p1) const {

    //coin flip
    bool coinFlip = (bool) std::rand() % 2;
    if (coinFlip){
        // pick a random index from the current point
        unsigned long randomCoOrdinate = std::rand() % keys.size();
        unsigned long index = keys[randomCoOrdinate];

        // Check if the index exists in the other point
        auto search1 = sparsePoint.find(index);
        auto search2 = p1.sparsePoint.find(index);
        if (search1 != p1.sparsePoint.end() )
            return (float) 0.5*std::abs(search1->second - search2->second);
        else
            return std::abs(search1->second);
    }
    else{
        // pick a random index from the current point
        unsigned long randomCoOrdinate = std::rand() % p1.keys.size();
        unsigned long index = p1.keys[randomCoOrdinate];

        // Check if the index exists in the other point
        auto search1 = sparsePoint.find(index);
        auto search2 = p1.sparsePoint.find(index);
        if (search2 != sparsePoint.end() )
            return (float) 0.5*std::abs(search1->second - search2->second);
        else
            return std::abs(search2->second);
    }

}



