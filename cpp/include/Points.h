//
// Created by Vivek Kumar Bagaria on 1/31/18.
//
#include <vector>
#include <iostream>
#include <memory>
#include <map>
#include <unordered_map>

#ifndef COMBINATORIAL_MAB_POINTS_H
#define COMBINATORIAL_MAB_POINTS_H

class BasePoint{
public:
    unsigned long vecSize;
    unsigned long dimensionPointer;

    BasePoint();

    virtual float distance(const BasePoint &p1) const;
    virtual std::pair<float, float> sampleDistance(const BasePoint &p1) const;
    unsigned long getVecSize() const;

};
class Point: public BasePoint {

public:
    std::vector <float> point;

    explicit Point();
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

    explicit SparsePoint();
    explicit SparsePoint(const std::unordered_map <unsigned long, float> &sp, unsigned long d);
};


class SparseL1Point : public SparsePoint {
public:
    explicit SparseL1Point(const std::unordered_map <unsigned long, float> &sp, unsigned long d);

    virtual float distance(const SparseL1Point &p1) const;
    virtual std::pair<float, float> sampleDistance(const SparseL1Point &p1, const unsigned sampleSize ) const;
};

class SquaredEuclideanPoint : public Point{
/* Points in Squared Euclidean space
 * */
public:
    explicit SquaredEuclideanPoint(const std::vector<float> &p);

    /*Computes the exact distance between two points.
     * Used only for debug purposes*/
    float distance(const SquaredEuclideanPoint &p1) const;
//    float distance(const std::shared_ptr<SquaredEuclideanPoint> p1) const;

    /*Picks a dimension of points pointed by the dimensionPointer and samples the distance
     * that dimension*/
    std::pair<float, float> sampleDistance(const SquaredEuclideanPoint &p1, const unsigned sampleSize)  const ;
};


class L1Point : public Point{
/* Points in Squared Euclidean space
 * */
public:
    explicit L1Point(const std::vector<float> &p);

    /*Computes the exact distance between two points.
     * Used only for debug purposes*/
    float distance(const L1Point &p1) const;

    /*Picks a dimension of points pointed by the dimensionPointer and samples the distance
     * that dimension*/
    std::pair<float, float> sampleDistance(const L1Point &p1, const unsigned sampleSize);
};

template <class templatePoint>
class GroupPoint: public BasePoint {
/*
 *
 * */
public:
    std::vector<std::shared_ptr<templatePoint> >  groupPoint;
    unsigned long noOfPoints;
    unsigned long groupID;

    explicit GroupPoint(){}
    explicit GroupPoint(const std::vector<std::shared_ptr<templatePoint> > &gp, unsigned d, unsigned long gid){
        groupPoint = gp;
        noOfPoints = gp.size();
        vecSize = d;
        groupID = gid;
    }

    float distance(const GroupPoint &gp1) const {
        float result(0);
        for(unsigned long i(0); i< noOfPoints; i++){
            for(unsigned long j(0); j< gp1.noOfPoints; j++){
//                std::cout << groupPoint.size() << std::endl;
                result += groupPoint[i]->distance((templatePoint)*gp1.groupPoint[j]);
            }
        }
        return result/(noOfPoints*gp1.noOfPoints);
    };

    std::pair<float, float> sampleDistance(const GroupPoint &gp1, const unsigned sampleSize) const {
        unsigned long i,j;
        i = std::rand() % noOfPoints;
        j = std::rand() % gp1.noOfPoints;
        std::pair<float, float> sample; ;
        sample = groupPoint[i]->sampleDistance((templatePoint)*gp1.groupPoint[j], sampleSize);
        return sample;
    };

    //For k-d trees
    inline unsigned int kdtree_get_point_count() const  {
        return noOfPoints ;
    };

    inline float kdtree_get_pt(const long idx, int dim) const {
        return groupPoint[idx].point[dim];
    };

    void addNewPoint(const std::shared_ptr<templatePoint> &np){
        groupPoint.push_back(np);
        noOfPoints ++;
    };
};



#endif //COMBINATORIAL_MAB_POINTS_H
