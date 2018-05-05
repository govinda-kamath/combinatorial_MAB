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
    std::string nameConfig = argv[1];
//    std::string nameConfig = "/Users/vivekkumarbagaria/Code/combinatorial_MAB/nominal.ini";
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
    int testnum = 500;
    int* indices = new int[testnum];
    for ( int i(0); i < testnum; i++) {
        indices[i] = std::rand() % n;
    }

    float** dist = new float*[testnum];
    for(int i = 0; i < testnum; ++i)
        dist[i] = new float[n];

    L2Oracle oracle(allPointsVec);

    std::cout << "Computing brute distances" << std::endl;
    for( int i(0); i<testnum; i++){
        int index = indices[i];
        for( int j(0); j<n; j++){
            float tmp = oracle(index,j);
            dist[i][j] = tmp;
        }
    }


    long totalCost = 0;
    std::cout << "Starting NNdescent for " << n << " points" << std::endl;
    std::cout << "S = " << S << std::endl;
    similarity::NNDescent<L2Oracle> nndescent(n ,k , S, oracle);

    std::cout << "Iterations " << std::endl;
    // Iterations
    for( int iter(1); iter< 30; iter++){
        std::chrono::system_clock::time_point timeStart = std::chrono::system_clock::now();
        nndescent.iterate(true);
        std::chrono::system_clock::time_point timeEnd = std::chrono::system_clock::now();
        totalCost += nndescent.getCost();
        long long int totalTime = std::chrono::duration_cast<std::chrono::milliseconds>
                (timeEnd - timeStart).count();


        std::cout << "Scanrate: " <<  totalCost/(n*n*0.5);

        float errorRate = 0;
        std::vector<similarity::KNN> nn = nndescent.getNN();
        for ( int i(0); i < testnum; i++){
            int item = indices[i];
            std::vector<int> v(dist[i], dist[i] + n);
            std::vector<int> sort_item = utils::sort_indexes(v);
            std::vector<int>  nndes;


            for( int j(0); j<k ; j++)
                nndes.push_back(nn[item][j].key);

            std::vector<int> v_intersection;

            std::set_union(sort_item.begin()+1, sort_item.begin()+k+1,
                           nndes.begin(), nndes.end(),  std::back_inserter(v_intersection));

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