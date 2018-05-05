//
// Created by Vivek Kumar Bagaria on 5/4/18.
//

#include "nndes.h"
#include "nndes-common.h"
#include "Arms.h"
#include "INIReader.h"
#include <vector>
#include <dlib/image_io.h>
#include <dlib/image_transforms.h>
#include "utils.h"

class L2Oracle {
public:
    std::vector<SquaredEuclideanPoint> pointsVec;
    float** dist; //memoization

    L2Oracle(std::vector<SquaredEuclideanPoint> pointsVec_) :
            pointsVec(pointsVec_) {
        long n = pointsVec.size();
        dist = new float*[n];
        for(long i = 0; i < n; ++i) {
            dist[i] = new float[n];
        }

        for(long i = 0; i < n; ++i) {
            for (long j = 0; j < n; ++j)
                dist[i][j] = -1;
        }
    }

    L2Oracle(std::vector<SquaredEuclideanPoint> pointsVec_, float** dist_) {
            pointsVec = pointsVec_;
            dist = dist_;
    }

    inline float operator()(unsigned long id1, unsigned long id2) const {
        if(dist[id1][id2] == -1){
            dist[id1][id2] = pointsVec[id1].distance(pointsVec[id2]);
            dist[id2][id1] = dist[id1][id2];
        }
        return dist[id1][id2];
    }
};


int main(int argc, char *argv[]) {
//    std::string nameConfig = argv[1];
    std::string nameConfig = "/Users/vivekkumarbagaria/Code/combinatorial_MAB/nominal.ini";
    INIReader reader(nameConfig);
    if (reader.ParseError() < 0) {
        std::cout << "Can't load "<< nameConfig << std::endl;
        return 1;
    }
    std::string directoryPath = reader.Get("path", "directory", "");
    std::string fileSuffix = reader.Get("path", "suffix", "");
    long n = (unsigned) reader.GetInteger("UCB", "n", -1);
    unsigned k = (unsigned) reader.GetInteger("UCB", "k", 50);
    unsigned S = (unsigned) reader.GetInteger("UCB", "S", 20);
    std::cout << "Running "<<k<< "-means for " << n<< " points" << std::endl;


    // Data
    std::vector<std::string>  pathsToImages;
    utils::getPathToFile(pathsToImages, directoryPath, fileSuffix);
    //Points
    std::vector<SquaredEuclideanPoint> allPointsVec;
    utils::vectorsToPoints(allPointsVec, pathsToImages, n);

    /*
     * Brute Method
     */
    float** dist = new float*[n];
    for(int i = 0; i < n; ++i)
        dist[i] = new float[n];

    L2Oracle oracle(allPointsVec);

//#define BruteComputed
#ifndef BruteComputed
    std::cout << "Computing brute distances" << std::endl;
    std::ofstream saveFile;
    saveFile.open("../NNdescentBrute"+std::to_string(n), std::ofstream::out | std::ofstream::trunc);
        std::vector<std::vector<int>> answers;
        for( int i(0); i<n; i++){
            dist[i][i] = 0;
            for( int j(0); j<i; j++){
                float tmp = oracle(i,j);
                dist[i][j] = tmp;
                dist[j][i] = tmp;
            }
        }
//        utils::serialize(saveFile, dist, n, n );
#else
        std::cout << "Reading off brute distances" << std::endl;

        std::ifstream saveFile;
        saveFile.open("NNdescentBrute"+std::to_string(n), std::ifstream::in);
        long rowtmp, coltmp;
        dist = utils::deserialize(saveFile, rowtmp, coltmp );
#endif


    long totalCost = 0;
    int testnum = 100;
    std::cout << "Starting NNdescent for " << n << " points" << std::endl;
    std::cout << "S = " << S << std::endl;
    similarity::NNDescent<L2Oracle> nndescent(n ,k , S, oracle);
//    std::cout << "Initial cost " << nndescent.getCost() << std::endl;
    for( int i(1); i< 30; i++){
        std::chrono::system_clock::time_point timeStart = std::chrono::system_clock::now();
        nndescent.iterate(true);
        std::chrono::system_clock::time_point timeEnd = std::chrono::system_clock::now();
        totalCost += nndescent.getCost();
        long long int totalTime = std::chrono::duration_cast<std::chrono::milliseconds>
                (timeEnd - timeStart).count();

//        std::cout << "Total cost " << totalCost/n << " with time " << totalTime <<  std::endl;
//        std::cout << "Estimated cost " << i*S*S*k*k << std::endl;
//        std::cout << "Ratio " <<  totalCost/(n*i*S*S*k*k+0.0) << std::endl;
        std::cout << "Scanrate: " <<  totalCost/(n*n*0.5);

        float errorRate = 0;
        std::vector<similarity::KNN> nn = nndescent.getNN();
        for ( int index(0); index < testnum; index++){
            int item = std::rand()%n;
            std::vector<int> v(dist[item], dist[item] + n);
            std::vector<int> sort_item = utils::sort_indexes(v);
            std::vector<int>  nndes;

//            std::cout << "Brute " << item <<" : ";
//            for( int j(1); j<k+1 ; j++){
//                std::cout << sort_item[j] << " ";
//            }
//            std::cout << std::endl;
//            std::cout << "NNdes " << item <<" : ";

            for( int j(0); j<k ; j++)
//            {
//                std::cout << nn[item][j].key << " ";
                nndes.push_back(nn[item][j].key);

//            }
//            std::cout << std::endl;
//            std::cout << "Common " << item <<" : ";

            std::vector<int> v_intersection;
            std::set_union(sort_item.begin()+1, sort_item.begin()+k+1,
                           nndes.begin(), nndes.end(),
                           std::back_inserter(v_intersection));

//            for(int n : v_intersection)
//                std::cout << n << ' ';
//            std::cout << std::endl;
//            std::cout << "Error: " << (v_intersection.size()-k)/(0.0+k) << std::endl;
            errorRate += (v_intersection.size()-k)/(0.0+k);
        }
        errorRate /= testnum;
        std::cout << "\t. Accuracy: " << 1- errorRate << "\n" << std::endl;
    }

    /*
    float fnnAvg, rnnAvg , fnnAvg2, rnnAvg2;
    fnnAvg = 0;
    rnnAvg = 0;
    std::vector<std::vector<int> > nn_new = nndescent.nn_new;
    std::vector<std::vector<int> > rnn_new = nndescent.rnn_new;

    for( int i(0); i<n ; i++){
        fnnAvg += nn_new[i].size();
        rnnAvg += rnn_new[i].size();
        fnnAvg2 += nn_new[i].size()*nn_new[i].size();
        rnnAvg2 += rnn_new[i].size()*rnn_new[i].size();
    }

    std::cout << "Average size of fnn is " << fnnAvg/n << std::endl;
    std::cout << "Sigma of fnn is " << fnnAvg2/n - (fnnAvg/n)*(fnnAvg/n) << std::endl;
    std::cout << "Average size of rnn is " << rnnAvg/n << std::endl;
    std::cout << "Sigma of rnn is " << rnnAvg2/n - (rnnAvg/n)*(rnnAvg/n) << std::endl;
    */


}