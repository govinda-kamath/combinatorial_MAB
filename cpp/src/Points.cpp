//
// Created by Vivek Kumar Bagaria on 1/31/18.
//

#include "Points.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <thread>
#include <dlib/image_io.h>


Point::Point(std::vector <float> p){
    point = p;
    vecSize = p.size();
}

unsigned long Point::getVecSize() const {
    return vecSize;
}
float Point::distance(const Point& p1) const {}
float Point::sampledDistance(const Point& p1) const {}

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

