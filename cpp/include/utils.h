//
// Created by Govinda Kamath on 2/1/18.
//

#include <vector>
#include <queue>
#include <string>

#ifndef COMBINATORIAL_MAB_UTILS_H
#define COMBINATORIAL_MAB_UTILS_H

#endif //COMBINATORIAL_MAB_UTILS_H

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
}