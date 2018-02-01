//
// Created by Vivek Kumar Bagaria on 1/31/18.
//
#include <vector>
#include <iostream>
#ifndef COMBINATORIAL_MAB_POINTS_H
#define COMBINATORIAL_MAB_POINTS_H


class Point {
/*A point is the base class that defines a point.
 * It has a vector that stores the representation of
 * the point in some high dimensional space.
 * This class also defines a distance function which
 * defines distances.
 * Points in different hilbert spaces are inherited
 * from this class.
 * */
public:
    std::vector <float> point;
    unsigned long vecSize;
    Point(std::vector <float> p);

    virtual float distance(const Point& p1) const;
    virtual float sampledDistance(const Point& p1) const;
    unsigned long getVecSize() const;
};

class SquaredEuclideanPoint : public Point{
/* Points in Squared Euclidean space
 * */
public:
    SquaredEuclideanPoint(std::vector<float> p);

    /*Computes the exact distance between two points.
     * Used only for debug purposes*/

    float distance(const SquaredEuclideanPoint& p1) const;

    /*Picks a dimension of points randomly and samples the distance
     * that dimension*/
    float sampledDistance(const SquaredEuclideanPoint& p1) const;
};



#endif //COMBINATORIAL_MAB_POINTS_H
