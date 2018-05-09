#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <random>
#include <thread>
#include <numeric>
#include <map>
#include <dlib/image_io.h>
#include <dlib/image_transforms.h>
#include "Points.h"
#include "INIReader.h"
#include "kmeans.h"


int main(int argc, char *argv[]){

    std::default_random_engine generator;
    std::normal_distribution<double> distribution1(0,1.0);
    std::normal_distribution<double> distribution2(0,2);
    std::cout << "Started!" << std::endl;
    std::string str;
    float tmp;

    long n = 20000;
    long m =10;
    std::vector<SquaredEuclideanPoint> dataMatrix;

    std::ifstream file("/Users/vivekkumarbagaria/Code/test_dataset/gene_data_MMI.txt", std::ios::binary);
    int i = 0;
    while (getline(file, str)) {
        std::stringstream ss(str);
        std::cout << i++ << " ";
        std::vector<float> tmpVec;
        for (long j = 0; j < m; j++){//5000
            ss >> tmp;
            tmpVec.push_back(tmp);
        }
        dataMatrix.push_back(SquaredEuclideanPoint(tmpVec));
        if(i>100)
            break;
    }
    std::cout << std::endl;

    std::vector<Arm2DMutualInformation<SquaredEuclideanPoint> > armsVec;
    for(unsigned i(1); i< m ; i ++) {
        for (unsigned j(0); j < m; j++) {
            std::vector<unsigned long> indices = {i, j};
            Arm2DMutualInformation<SquaredEuclideanPoint> arm(1, dataMatrix, indices);
            armsVec.push_back(arm);
        }
    }

    //Arm1
    std::vector<SquaredEuclideanPoint> arm1vec, arm2vec;
    for (int i=0; i<n; ++i) {

        float number1 = (float) distribution1(generator);
        float number2 = (float) distribution1(generator);
        std::vector<float> tmp1 = {number1, number1+number2};
        arm1vec.push_back(SquaredEuclideanPoint(tmp1));


        number1 = (float) distribution2(generator);
        number2 = (float) distribution2(generator);
        std::vector<float> tmp2 = {number1, 2*number1+number2};
        arm2vec.push_back(SquaredEuclideanPoint(tmp2));
    }
    std::ofstream file2;
    file2.open("/Users/vivekkumarbagaria/Code/combinatorial_MAB/gaussians", std::ofstream::out | std::ofstream::trunc);

    for (long i = 0; i < n; i++) {
        file2 << arm1vec[i].point[0] << "," << arm1vec[i].point[1] <<std::endl;
//        return arr;
    }

    Arm2DMutualInformation<SquaredEuclideanPoint> arm1(1, arm1vec, indices);
    Arm2DMutualInformation<SquaredEuclideanPoint> arm2(1, arm2vec, indices);

    std::vector<Arm2DMutualInformation<SquaredEuclideanPoint>> allArmsVec;
    allArmsVec.push_back(arm1);
    allArmsVec.push_back(arm2);

    for(int i=0;i < 5; i++){
        arm1.pullArm(0, 0, 1, true, (int) n/4, -1);
        arm2.pullArm(0, 0, 1, true, (int) n/4, -1);
        std::cout << "Arm 1 " << arm1.lowerConfidenceBound << " "
                  << arm1.estimateOfMean << " " << arm1.upperConfidenceBound << " " << "0.346" << std::endl;
        std::cout << "Arm 2 " << arm2.lowerConfidenceBound << " "
                  << arm2.estimateOfMean << " " << arm2.upperConfidenceBound << " " << 0.80 << std::endl;
        std::cout << arm2.estimateOfMean - arm1.estimateOfMean << " 0.454\n" <<std::endl;

    }

}
