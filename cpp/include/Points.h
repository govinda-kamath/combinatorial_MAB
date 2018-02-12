//
// Created by Vivek Kumar Bagaria on 1/31/18.
//
#include <vector>
#include <iostream>
#include <unordered_map>

#ifndef COMBINATORIAL_MAB_POINTS_H
#define COMBINATORIAL_MAB_POINTS_H

class BasePoint{
public:
    unsigned long vecSize;

    virtual float distance(const BasePoint &p1) const;
    virtual float sampledDistance(const BasePoint &p1) const;
    unsigned long getVecSize() const;

};
class Point: public BasePoint {

public:
    std::vector <float> point;
    explicit Point(const std::vector <float> &p);
};


class SparsePoint: public BasePoint {
/*A point is the base class that defines a point.
 * It has a vector that stores the representation of
 * the point in some high dimensional space.
 * This class also defines a distance function which
 * defines distances.
 * Points in different hilbert spaces are inherited
 * from this class.
 * */
public:
    std::unordered_map <unsigned long, float>  sparsePoint;
    std::vector<unsigned long> keys;

    explicit SparsePoint(const std::unordered_map <unsigned long, float> &sp, unsigned long d);
};


class SparseL1Point : public SparsePoint {
public:
    explicit SparseL1Point(const std::unordered_map <unsigned long, float> &sp, unsigned long d);

    virtual float distance(const SparseL1Point &p1) const;
    virtual float sampledDistance(const SparseL1Point &p1) const;
};

class SquaredEuclideanPoint : public Point{
/* Points in Squared Euclidean space
 * */
public:
    explicit SquaredEuclideanPoint(const std::vector<float> &p);

    /*Computes the exact distance between two points.
     * Used only for debug purposes*/
    float distance(const SquaredEuclideanPoint &p1) const;

    /*Picks a dimension of points randomly and samples the distance
     * that dimension*/
    float sampledDistance(const SquaredEuclideanPoint &p1) const;
};


class L1Point : public Point{
/* Points in Squared Euclidean space
 * */
public:
    explicit L1Point(const std::vector<float> &p);

    /*Computes the exact distance between two points.
     * Used only for debug purposes*/
    float distance(const L1Point &p1) const;

    /*Picks a dimension of points randomly and samples the distance
     * that dimension*/
    float sampledDistance(const L1Point &p1) const;
};

//template <class templatePoint>
//class GroupPoint: public BasePoint {
///*A point is the base class that defines a point.
// * It has a vector that stores the representation of
// * the point in some high dimensional space.
// * This class also defines a distance function which
// * defines distances.
// * Points in different hilbert spaces are inherited
// * from this class.
// * */
//public:
//    std::vector<templatePoint>  groupPoint;
//    std::vector<unsigned long> noOfPoints;
//
//    explicit GroupPoint(std::unordered_map <unsigned long, float> sp, unsigned long d);
//};


#endif //COMBINATORIAL_MAB_POINTS_H
