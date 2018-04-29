//
// Created by Govinda Kamath on 2/5/18.
//

#ifndef COMBINATORIAL_MAB_KMEANS_H
#define COMBINATORIAL_MAB_KMEANS_H

#endif //COMBINATORIAL_MAB_KMEANS_H
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <thread>
#include <numeric>
#include <map>
#include <dlib/image_io.h>
#include <dlib/image_transforms.h>
#include "Points.h"
#include "Arms.h"
#include "../deprecated/UCB.h"
#include <stdexcept>
#include "utils.h"
#include "Knn.h"

template <class templatePoint>
class Kmeans : public Knn<templatePoint>{
public:
    Kmeans( std::vector<templatePoint>& pVecL, std::vector<templatePoint>& pVecR,
            unsigned noOfInitialPulls, float deltaAccuracy ):
            Knn<templatePoint>(pVecL,  pVecR, 1,  noOfInitialPulls,  deltaAccuracy ){}

    using Knn<templatePoint>::run;

    using Knn<templatePoint>::saveAnswers;

};