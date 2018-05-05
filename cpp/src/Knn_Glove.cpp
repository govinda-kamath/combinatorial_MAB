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
#include "INIReader.h"
#include "kmeans.h"


int main(int argc, char *argv[]){

    std::string nameConfig = argv[1];
    long startIndex(atol(argv[2])); // Start index
    long endIndex(atol(argv[3])); // End index
    INIReader reader(nameConfig);
    if (reader.ParseError() < 0) {
        std::cout << "Can't load "<< nameConfig << std::endl;
        return 1;
    }
    unsigned numberOfInitialPulls = (unsigned) reader.GetInteger("UCB", "numberOfInitialPulls_knn", 100);
    std::string saveFileFolder =reader.Get("path", "saveFileFolderKnnGlove", "experiments/knn/tmp");
    std::string directoryPath = reader.Get("path", "directoryGlove", "/Users/vivekkumarbagaria/Code/test_dataset/glove.840B.300d.txt");

    unsigned k = (unsigned) reader.GetInteger("UCB", "k", 50);
    float delta = (float) reader.GetReal("UCB", "delta", 0.1);
    unsigned sampleSize = (unsigned) reader.GetInteger("UCB", "sampleSize", 32);
    long n = (unsigned) reader.GetInteger("UCB", "n", -1);
    delta = delta / (k);

    std::ifstream file(directoryPath);


    std::vector<std::vector<float>> data;
    std::vector<std::string> name;
    for(unsigned i(0); i<n; i++){
        std::vector<float> tmpV(300);
        data.push_back(tmpV);
        name.push_back("");
    }
    std::string tmp;

    for(int row = 0; row < n; ++row)
    {
        std::string line;
        std::getline(file, line);
        std::stringstream iss(line);

        for (int col = 0; col < 301; ++col)
        {
            std::string val;
            std::getline(iss, val, ' ');

            std::stringstream convertor(val);
            if (col==0){
                convertor >> name[row];
                continue;
            }
            convertor >> data[row][col];
//            std::cout << row << " " << col << " " << data[row][col] << std::endl;
        }
    }
    std::vector<SquaredEuclideanPoint> pointsVec;
    utils::vectorsToPoints(pointsVec, data);


    std::chrono::system_clock::time_point loopTimeStart = std::chrono::system_clock::now();
    Knn<SquaredEuclideanPoint> knn(pointsVec, k, numberOfInitialPulls, delta, sampleSize, saveFileFolder);
    std::vector<unsigned long> indices(endIndex-startIndex);
    std::iota(indices.begin(), indices.end(), startIndex);
    std::cout << "Running" << std::endl;
    knn.run(indices);
    std::chrono::system_clock::time_point loopTimeEnd = std::chrono::system_clock::now();

    std::cout << "Average time (ms) "
              << std::chrono::duration_cast<std::chrono::milliseconds>(loopTimeEnd - loopTimeStart).count()/
                 (endIndex-startIndex) << std::endl;




}
