//
// Created by Govinda Kamath on 2/1/18.
//

#include <vector>
#include <algorithm>
#include <queue>
#include <string>
#include <unordered_map>
#include <utility>

#ifndef COMBINATORIAL_MAB_UTILS_H
#define COMBINATORIAL_MAB_UTILS_H


namespace utils{
    //Reads the protected container object of a priority queue to give access to its elements
    //A beautiful hack picked up from https://stackoverflow.com/a/12886393/3377314
    template <class T, class S, class C>
    S& Container(std::priority_queue<T, S, C>& q) {
        struct HackedQueue : private std::priority_queue<T, S, C> {
            static S& Container(std::priority_queue<T, S, C>& q) {
                return q.*&HackedQueue::c;
            }
        };
        return HackedQueue::Container(q);
    }

    void readImageAsVector (std::string filePath, std::vector<float> &imageVec);

    void getPathToFile(std::vector<std::string> & pathsToImages, const std::string directoryPath,
    const std::string  fileSuffix);

    template <class templatePoint>
    void vectorsToPoints(std::vector<templatePoint> &pointsVec,
                         std::vector<std::string>  &pathsToImages){
        for  (unsigned long i(0); i < pathsToImages.size(); i++) {
            std::vector<float> tmpVec;
            readImageAsVector(pathsToImages[i],tmpVec);
            templatePoint tmpPoint(tmpVec);
            pointsVec.push_back(tmpPoint);
            if (i%10000 == 9999){
                std::cout << i+1 << " points read." << std::endl;
            }
        }
    }

    template <class templatePoint, class dataType>
    void vectorsToPoints(std::vector<templatePoint> &pointsVec, std::vector<std::vector<dataType> > &dataMatrix) {
        //Each row in the data-matrix is a point

        for (unsigned long i(0); i < dataMatrix.size(); i++) {

            std::vector<float> tmpVec;
            for (unsigned long j(0); j < dataMatrix[0].size(); j++) {
                tmpVec.push_back((float)dataMatrix[i][j]);
            }

            templatePoint tmpPoint(tmpVec);
            pointsVec.push_back(tmpPoint);

        }
    }

    template <class templatePoint>
    void unorderedMapToPoints(std::vector<templatePoint> &pointsVec,
                         std::vector<std::unordered_map<unsigned long, float> > &sparseDataMatrix,
                              unsigned long dimension) {

        //Each row in the data-matrix is a point

        for (unsigned long i(0); i < sparseDataMatrix.size(); i++) {
            templatePoint tmpPoint(sparseDataMatrix[i], dimension);
            pointsVec.push_back(tmpPoint);
        }
    }

    struct ArmConditions{
        unsigned long numberOfPulls;
        float sumOfPulls;
        float sumOfSquaresOfPulls;
        float trueMeanValue;

        ArmConditions () : numberOfPulls(0), sumOfPulls(0.0), sumOfSquaresOfPulls(0.0), trueMeanValue(INFINITY){}
        ArmConditions (unsigned long nP, float sumPulls, float sumSquare, float tMeanValue) : numberOfPulls(nP),
                       sumOfPulls(sumPulls), sumOfSquaresOfPulls(sumSquare), trueMeanValue(tMeanValue){}
        void update (unsigned long nP, float sumPulls, float sumSquare){
            numberOfPulls = nP;
            sumOfPulls = sumPulls;
            sumOfSquaresOfPulls = sumSquare;
        }
    };

    // Copied it from
    // https://stackoverflow.com/a/32685618/803072
    struct pair_hash {
        template <class T1, class T2>
        std::size_t operator () (const std::pair<T1,T2> &p) const {
            auto h1 = std::hash<T1>{}(p.first);
            auto h2 = std::hash<T2>{}(p.second);

            // Mainly for demonstration purposes, i.e. works but is overly simple
            // In the real world, use sth. like boost.hash_combine
            return h1 ^ h2;
        }
    };
}
#endif //COMBINATORIAL_MAB_UTILS_H
