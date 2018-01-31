#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <thread>
#include <dlib/image_io.h>
#include <dlib/image_transforms.h>
#include <glob.h>

#define DEBUG

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
    explicit Point(std::vector <float> p){
        point = p;
    }
    ~Point () {
//        delete &point;
    }

    virtual float distance(const Point& p1){}
    virtual float sampledDistance(const Point& p1){}
};

class SquaredEuclideanPoint : public Point{
/* Points in Squared Euclidean space
 * */
public:
    explicit SquaredEuclideanPoint(std::vector<float> p) : Point(p){}

    /*Computes the exact distance between two points.
     * Used only for debug purposes*/
    float distance(const SquaredEuclideanPoint& p1){

        assert(("Sizes do not match", point.size() == p1.point.size()));

        float result(0);

        std::vector<float>::const_iterator pIt = point.begin();
        std::vector<float>::const_iterator p1It = p1.point.begin();
        for (; p1It != p1.point.end() && pIt  != point.end(); ++p1It, ++pIt){
            result += (*p1It-*pIt)*(*p1It-*pIt);
        }
        return result;
    }

    /*Picks a dimension of points randomly and samples the distance
     * that dimension*/
    float sampledDistance(const SquaredEuclideanPoint& p1){
        assert(("Sizes do not match", point.size() == p1.point.size()));

        unsigned vecSize = p1.point.size();
        float result;
        unsigned randomCoOrdinate;
        randomCoOrdinate = std::rand() % vecSize;
        return (point[randomCoOrdinate] - p1.point[randomCoOrdinate])
               *(point[randomCoOrdinate] - p1.point[randomCoOrdinate]);
    }
};



template <class templatePoint>
class Arm {
    /* A base class of arms */
public:
    unsigned long numberOfPulls;
    float sumOfPulls;
    float upperConfidenceBound;
    float lowerConfidenceBound;
    float estimateOfMean;
    float estimateOfSecondMoment;
    float SumOfSquaresOfPulls;
    unsigned long dimension;
    templatePoint * point;
    unsigned long id;

    Arm(){
        numberOfPulls = 0;
        sumOfPulls = 0.0;
        upperConfidenceBound = INFINITY;
        lowerConfidenceBound = -INFINITY;
        SumOfSquaresOfPulls = 0.0;
        estimateOfMean = NAN;
        estimateOfSecondMoment = NAN;
    }

    Arm(unsigned long armNumber, templatePoint p) : Arm() {
        id = armNumber;
        point = new templatePoint(p.point);
        dimension = point->point.size();
    }

    ~Arm(){
//        delete point; //no clue why this errors out
    }

    void printArm(){
        std::cout << "Number of pulls" << numberOfPulls
                                       << std::endl;
    }

    friend bool operator> (const Arm& l, const Arm& r)
        {
            return l.lowerConfidenceBound > r.lowerConfidenceBound;
        }

    float updateConfidenceIntervals(float globalSigma, float logDeltaInverse){

        float localSigma, intervalWidth;
        localSigma = globalSigma; //Todo: update sigma to new local value
        intervalWidth = std::sqrt((localSigma * localSigma * logDeltaInverse)/numberOfPulls);
        upperConfidenceBound = estimateOfMean + intervalWidth;
        lowerConfidenceBound = std::max((float)0.0, estimateOfMean - intervalWidth);
    }

    float pullArm(templatePoint &p1, float globalSigma,
                            float logDeltaInverse, bool update = true) {
        float sample;

        if (numberOfPulls >= dimension){
            sample = trueMean(p1);
            numberOfPulls += dimension;
            estimateOfMean = sample;
            upperConfidenceBound = estimateOfMean;
            lowerConfidenceBound = estimateOfMean;
            sumOfPulls += sample*dimension;
            SumOfSquaresOfPulls += std::pow(sample,2)*dimension;
            estimateOfSecondMoment = sample*sample;
        }
        else {
            sample = point->sampledDistance(p1);

            numberOfPulls++;
            sumOfPulls += sample;
            estimateOfMean = sumOfPulls / numberOfPulls;
            SumOfSquaresOfPulls += sample * sample;
            estimateOfSecondMoment = SumOfSquaresOfPulls / numberOfPulls;

            if (update)
                updateConfidenceIntervals(globalSigma, logDeltaInverse);
        }

        return sample;
    }

    float trueMean(const templatePoint &p1){
        return point->distance(p1)/p1.point.size();
    }

};


template <class templatePoint>
class ArmKNN : public Arm<templatePoint>{
    /*
     * Arms for k-Nearest Neighbours*/
public:
    templatePoint *fixedPoint;

    ArmKNN(unsigned long id, templatePoint p) : Arm<templatePoint>(id, p) {}

    ArmKNN(unsigned long id, templatePoint p, templatePoint fixPoint) : Arm<templatePoint>(id, p) {
        fixedPoint = new templatePoint(fixPoint.point);
    }

    using Arm<templatePoint>::pullArm;
    float pullArm(float globalSigma, float logDeltaInverse, bool update = true){
        return pullArm(*fixedPoint, globalSigma, logDeltaInverse, update);
    }

    using Arm<templatePoint>::trueMean;
    float trueMean(){
        return trueMean(*fixedPoint);
    }
};


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


template <class templateArm>
class UCB{
    /* UCB for the general case*/
public:
    unsigned numberOfArms;
    std::vector<templateArm> armsContainer;
    std::priority_queue<templateArm, std::vector<templateArm>, std::greater<templateArm> > arms;
    float logDeltaInverse;

    float globalSigma;
    float globalNumberOfPulls;
    float globalSumOfPulls;
    float globalSumOfSquaresOfPulls;
    std::mutex initializeMutex;
    UCB(std::vector<templateArm> &armsVec, float delta){

        armsContainer = armsVec;
        numberOfArms = armsContainer.size();

        logDeltaInverse = std::log(1/delta);

        globalSigma = NAN;
        globalNumberOfPulls = 0;
        globalSumOfPulls = 0;
        globalSumOfSquaresOfPulls = 0;

    }

    void initialiseFewArm(unsigned long armIndexStart, unsigned long armIndexEnd, int numberOfInitialPulls){

        float localSumOfPulls, localSumOfSquaresOfPulls;
        localSumOfPulls = 0;
        localSumOfSquaresOfPulls = 0;
        // Pulling an arm numberOfInitialPulls times
        for (unsigned long armIndex = armIndexStart; armIndex< armIndexEnd; armIndex++) {
            for (unsigned i = 0; i < numberOfInitialPulls; i++) {
                float observedSample(0);
                observedSample = armsContainer[armIndex].pullArm(0, 0, false);
                localSumOfPulls += observedSample;
                localSumOfSquaresOfPulls += observedSample * observedSample;

            }
        }
        // locking the global variables which have to be updated
        {
            std::lock_guard<std::mutex> guard(initializeMutex);
            globalSumOfPulls += localSumOfPulls;
            globalSumOfSquaresOfPulls += localSumOfSquaresOfPulls;
        }
    }
    void initialise(int numberOfInitialPulls = 100){

        unsigned numberOfThreads = 1 ;//std::thread::hardware_concurrency();
        std::cout << numberOfThreads << " number of threads used in each batch of initialization.\n";
        std::vector<std::thread> initThreads(numberOfThreads);

        unsigned long chunkSize = (numberOfArms/numberOfThreads);
        for(unsigned t = 0; t < numberOfThreads; t++){
            unsigned long armIndexStart = t*chunkSize;
            unsigned long armIndexEnd = (t+1)*chunkSize;

//            std::cout << "Starting thread for arm " << armIndexStart
//                      << " to " << armIndexEnd << std::endl;
            initThreads[t] = std::thread(&UCB::initialiseFewArm, this, armIndexStart, armIndexEnd, numberOfInitialPulls);
        }

        for(unsigned t = 0; t < numberOfThreads; t++){
//            std::cout << "Joining thread for group " << t << std::endl;
            initThreads[t].join();
        }




        globalNumberOfPulls = numberOfInitialPulls*numberOfArms;
        globalSigma = std::sqrt((globalSumOfSquaresOfPulls/globalNumberOfPulls -
                                    std::pow(globalSumOfPulls/globalNumberOfPulls,2)));
#ifdef DEBUG
        std::cout << "Sigma after initialization " << globalSigma <<std::endl;
        std::cout << "Mean after initialization "
                  << globalSumOfPulls/globalNumberOfPulls <<std::endl;
        std::cout << "Mean Squared after initialization "
                  << std::pow(globalSumOfPulls/globalNumberOfPulls,2) <<std::endl;
        std::cout << "Second Moment after initialization "
                  << globalSumOfSquaresOfPulls/globalNumberOfPulls <<std::endl;
#endif

        for (unsigned long index = 0; index < numberOfArms; index++){

            armsContainer[index].updateConfidenceIntervals(globalSigma, logDeltaInverse);
            arms.push(armsContainer[index]);
        }

    }

    void runUCB(unsigned long maxIterations){
        for (unsigned long i(0); i < maxIterations; i++){
            bool stop;
            stop = iterationOfUCB();
            if (stop){
                break;
            }
        }
    }

    bool iterationOfUCB(){
        /*An iteration of UCB*/
        templateArm bestArm = arms.top();
        arms.pop();
        templateArm secondBestArm = arms.top();
        float UCBofBestArm, LCBofBestArm, LCBofSecondBestArm;
        UCBofBestArm = bestArm.upperConfidenceBound;
        LCBofBestArm = bestArm.lowerConfidenceBound;
        LCBofSecondBestArm = secondBestArm.lowerConfidenceBound;
        if (UCBofBestArm < LCBofSecondBestArm){
            //Checking if UCB should stop
            arms.push(bestArm);
#ifdef DEBUG
            std::cout << "stopping UCB "<< UCBofBestArm << std::endl;
            std::cout << "stopping LCB "<< LCBofSecondBestArm << std::endl;
#endif
            return true;
        } else {
            float sample;
            sample = bestArm.pullArm(globalSigma, logDeltaInverse);
            if (UCBofBestArm == LCBofBestArm){
                //Checking if the best arm is being computed in full and updating
                //things accordingly.
                unsigned long dimension(bestArm.dimension);
                globalNumberOfPulls += dimension;
                globalSumOfPulls += sample*dimension;
                globalSumOfSquaresOfPulls += std::pow(sample, 2)*dimension;
            }
            else {
                globalSumOfPulls += sample;
                globalSumOfSquaresOfPulls += std::pow(sample, 2);
                globalNumberOfPulls++;
            }
            globalSigma = std::sqrt((globalSumOfSquaresOfPulls / globalNumberOfPulls -
                                     std::pow(globalSumOfPulls / globalNumberOfPulls, 2)));
            arms.push(bestArm);
        }
        return  false;
    }
};

void readImageAsVector (std::string filePath, std::vector<float> &imageVec) {

    dlib::array2d <dlib::rgb_pixel> imageRGB;
    dlib::load_image(imageRGB, filePath.c_str());
    unsigned  numColumns(imageRGB.nc()), numRows(imageRGB.nr());
    unsigned numPixels(numColumns*numRows);
    unsigned vecLength(numPixels*3);

//    std::cout << numColumns <<"\t" << numRows <<
//              "\t" << numPixels << "\t" << vecLength << std::endl;

    if (imageVec.size() != vecLength){

//        std::cout << "initialising" << std::endl;
        imageVec.clear();
        imageVec.reserve(vecLength);

        for (unsigned i(0); i < numRows; i++){
            for (unsigned j(0); j < numColumns; j++) {
                imageVec.push_back((float) imageRGB[i][j].red);
            }
        }



        for (unsigned i(0); i < numRows; i++){
            for (unsigned j(0); j < numColumns; j++) {
                imageVec.push_back((float) imageRGB[i][j].blue);
            }
        }

        for (unsigned i(0); i < numRows; i++){
            for (unsigned j(0); j < numColumns; j++) {
                imageVec.push_back((float) imageRGB[i][j].green);
            }
        }

    } else {

        for (unsigned i(0); i < numRows; i++) {
            for (unsigned j(0); j < numColumns; j++) {
                imageVec[i * numRows + j] = (float) imageRGB[i][j].red;
            }
        }

        for (unsigned i(0); i < numRows; i++) {
            for (unsigned j(0); j < numColumns; j++) {
                imageVec[numPixels + i * numRows + j] = (float) imageRGB[i][j].blue;
            }
        }

        for (unsigned i(0); i < numRows; i++) {
            for (unsigned j(0); j < numColumns; j++) {
                imageVec[2 * numPixels + i * numRows + j] = (float) imageRGB[i][j].green;
            }
        }
    }
}


int main(int argc, char *argv[]){

    std::cout << "We have entered " << argc
         << " arguments." << std::endl;

    std::vector<SquaredEuclideanPoint> pointsVec;
    std::vector<ArmKNN<SquaredEuclideanPoint> > armsVec;
    int numberOfInitialPulls(200);
    float delta(0.001);

    if(argc != 4) {

//        std::string filePathTest("/Users/govinda/Code/combinatorial_MAB/val_27.JPEG");
//        readImageAsVector(filePathTest,tmpVectorImage);
//        std::cout << "length of read vector " << tmpVectorImage.size() << std::endl;
//        readImageAsVector(filePathTest,tmpVectorImage);
//        std::cout << "length of read vector " << tmpVectorImage.size() << std::endl;

        glob_t glob_result;
        std::vector<float> tmpVec;

        std::string directoryPath("/Users/govinda/Code/combinatorial_MAB/test_dataset/tiny-imagenet-200/val/images");
        std::string filePrefix("/val_");
        std::string fileSuffix(".JPEG");
        unsigned long fileNumber(0);
        std::string search_name;

        std::vector<std::string> pathsToImages;
        search_name = directoryPath + filePrefix + std::to_string(fileNumber) + fileSuffix;
        std::cout << search_name << std::endl;

        glob(search_name.c_str(),GLOB_TILDE,NULL,&glob_result);

        clock_t timeRead = clock();
        while (glob_result.gl_pathc != 0){
//            std::cout << std::string(glob_result.gl_pathv[0]) << std::endl;

            pathsToImages.push_back(std::string(glob_result.gl_pathv[0]));
            fileNumber ++;
            search_name = directoryPath + filePrefix + std::to_string(fileNumber) + fileSuffix;
            glob(search_name.c_str(),GLOB_TILDE,NULL,&glob_result);
//            std::cout << "Number of files " << glob_result.gl_pathc << std::endl;
        }

        unsigned long pointIndex(0);

        for  (unsigned long i(0); i < pathsToImages.size(); i++) {
            float tmpValue;

            std::vector<float> tmpVec;
            readImageAsVector(pathsToImages[i],tmpVec);

            SquaredEuclideanPoint tmpPoint(tmpVec);
            pointsVec.push_back(tmpPoint);
            pointIndex++;

            if (pointIndex%1000 == 999){
                std::cout << pointIndex+1 << " points read." << std::endl;
            }
        }

        for (unsigned i(1); i < pointsVec.size(); i++) {
            ArmKNN<SquaredEuclideanPoint> tmpArm(i - 1, pointsVec[i], pointsVec[0]);
            armsVec.push_back(tmpArm);
        }

        UCB<ArmKNN<SquaredEuclideanPoint> > UCB1(armsVec, delta);
        std::cout << "Reading time (ms)" << 1000 * (clock() - timeRead) / CLOCKS_PER_SEC << std::endl;

        clock_t timeInitialize = clock();
        UCB1.initialise(numberOfInitialPulls);
        std::cout << "Initializing time (ms)" << 1000 * (clock() - timeInitialize) / CLOCKS_PER_SEC << std::endl;

        std::vector<ArmKNN<SquaredEuclideanPoint> > &hackedArmsVec = Container(UCB1.arms);
#ifdef DEBUG
//        for (unsigned i = 0; i < 10; i++) {
//            std::cout << hackedArmsVec[i].id << " True Mean= " << armsVec[hackedArmsVec[i].id].trueMean()
//                      << " sigma = "
//                      << std::sqrt((hackedArmsVec[i].SumOfSquaresOfPulls / hackedArmsVec[i].numberOfPulls -
//                                    std::pow(hackedArmsVec[i].sumOfPulls / hackedArmsVec[i].numberOfPulls, 2)))
//                      << " estimate = " << hackedArmsVec[i].estimateOfMean << " total pulls="
//                      << hackedArmsVec[i].numberOfPulls << std::endl;
//        }

        std::cout << "average pull " << UCB1.globalNumberOfPulls / armsVec.size() << std::endl;
        std::cout << "sigma " << UCB1.globalSigma << std::endl;
        std::cout << "best arm's estimate " << UCB1.arms.top().estimateOfMean << std::endl;
        std::cout << UCB1.arms.top().id << std::endl;
#endif
                clock_t timeIterate = clock();
        UCB1.runUCB(10000000000);
        std::cout << "Iteration time (ms) " << 1000 * (clock() - timeIterate) / CLOCKS_PER_SEC << std::endl;

//    std::vector<ArmKNN<SquaredEuclideanPoint> > &hackedArmsVec = Container(UCB1.arms);

        clock_t timeTrueMean = clock();

        for (unsigned i = 0; i < hackedArmsVec.size(); i++) {
            float tmp = armsVec[hackedArmsVec[i].id].trueMean();
        }
        std::cout << "True Mean time (ms) " << 1000 * (clock() - timeTrueMean) / CLOCKS_PER_SEC << std::endl;

#ifdef DEBUG
        for (unsigned i = 0; i < hackedArmsVec.size(); i++) {
            std::cout << hackedArmsVec[i].id << " True Mean= " << armsVec[hackedArmsVec[i].id].trueMean()
                      << " sigma = "
                      << std::sqrt((hackedArmsVec[i].SumOfSquaresOfPulls / hackedArmsVec[i].numberOfPulls -
                                    std::pow(hackedArmsVec[i].sumOfPulls / hackedArmsVec[i].numberOfPulls, 2)))
                      << " estimate = " << hackedArmsVec[i].estimateOfMean
                      << " lcb = " << hackedArmsVec[i].lowerConfidenceBound
                      << " ucb = " << hackedArmsVec[i].upperConfidenceBound
                      << " total pulls="
                      << hackedArmsVec[i].numberOfPulls << std::endl;
        }

        std::cout << "average pull " << UCB1.globalNumberOfPulls / armsVec.size() << std::endl;
        std::cout << "sigma " << UCB1.globalSigma << std::endl;
        std::cout << "best arm's estimate " << UCB1.arms.top().estimateOfMean << std::endl;
        std::cout << UCB1.arms.top().id << std::endl;
#endif



    } else {
//        std::string filePath(argv[1]), line;
//        int numberOfInitialPulls(atoi(argv[2]));
//        float delta(atof(argv[3]));
//
////    filePath = "/Users/govinda/Code/combinatorial_MAB/test_dataset/1000_images.txt";
////    filePath ="/data/MAB/work/dataset/test_dataset/basic_io_dataset/10k_images.txt";
//        std::fstream fileReader(filePath.c_str());
//        unsigned long pointIndex(0);
//
//        std::vector<SquaredEuclideanPoint> pointsVec;
//        std::vector<ArmKNN<SquaredEuclideanPoint> > armsVec;
//
//        clock_t timeRead = clock();
//
//        while (std::getline(fileReader, line)) {
//            float tmpValue;
//
//            std::vector<float> tmpVec;
//            std::stringstream ss(line);
//            while (ss >> tmpValue) {
//                tmpVec.push_back(tmpValue);
//            }
//            SquaredEuclideanPoint tmpPoint(tmpVec);
//
//            pointsVec.push_back(tmpPoint);
//            pointIndex++;
//
//        }
//
//
//        for (unsigned i(1); i < pointsVec.size(); i++) {
//            ArmKNN<SquaredEuclideanPoint> tmpArm(i - 1, pointsVec[i], pointsVec[0]);
//            armsVec.push_back(tmpArm);
//        }
//
//
//        UCB<ArmKNN<SquaredEuclideanPoint> > UCB1(armsVec, delta);
//        std::cout << "Reading time (ms)" << 1000 * (clock() - timeRead) / CLOCKS_PER_SEC << std::endl;
//        clock_t timeInitialize = clock();
//        UCB1.initialise(numberOfInitialPulls);
//        std::cout << "Initializing time (ms)" << 1000 * (clock() - timeInitialize) / CLOCKS_PER_SEC << std::endl;
//
//        std::vector<ArmKNN<SquaredEuclideanPoint> > &hackedArmsVec = Container(UCB1.arms);
//#ifdef DEBUG
//        for (unsigned i = 0; i < hackedArmsVec.size(); i++) {
//            std::cout << hackedArmsVec[i].id << " True Mean= " << armsVec[hackedArmsVec[i].id].trueMean()
//                      << " sigma = "
//                      << std::sqrt((hackedArmsVec[i].SumOfSquaresOfPulls / hackedArmsVec[i].numberOfPulls -
//                                    std::pow(hackedArmsVec[i].sumOfPulls / hackedArmsVec[i].numberOfPulls, 2)))
//                      << " estimate = " << hackedArmsVec[i].estimateOfMean << " total pulls="
//                      << hackedArmsVec[i].numberOfPulls << std::endl;
//        }
//
//        std::cout << "average pull " << UCB1.globalNumberOfPulls / armsVec.size() << std::endl;
//        std::cout << "sigma " << UCB1.globalSigma << std::endl;
//        std::cout << "best arm's estimate " << UCB1.arms.top().estimateOfMean << std::endl;
//        std::cout << UCB1.arms.top().id << std::endl;
//#endif
//        clock_t timeIterate = clock();
//        UCB1.runUCB(10000000000);
//        std::cout << "Iteration time (ms) " << 1000 * (clock() - timeIterate) / CLOCKS_PER_SEC << std::endl;
//
////    std::vector<ArmKNN<SquaredEuclideanPoint> > &hackedArmsVec = Container(UCB1.arms);
//        clock_t timeTrueMean = clock();
//
//        for (unsigned i = 0; i < hackedArmsVec.size(); i++) {
//            float tmp = armsVec[hackedArmsVec[i].id].trueMean();
//        }
//        std::cout << "True Mean time (ms) " << 1000 * (clock() - timeTrueMean) / CLOCKS_PER_SEC << std::endl;
//
//#ifdef DEBUG
//        for (unsigned i = 0; i < hackedArmsVec.size(); i++) {
//            std::cout << hackedArmsVec[i].id << " True Mean= " << armsVec[hackedArmsVec[i].id].trueMean()
//                      << " sigma = "
//                      << std::sqrt((hackedArmsVec[i].SumOfSquaresOfPulls / hackedArmsVec[i].numberOfPulls -
//                                    std::pow(hackedArmsVec[i].sumOfPulls / hackedArmsVec[i].numberOfPulls, 2)))
//                      << " estimate = " << hackedArmsVec[i].estimateOfMean
//                      << " lcb = " << hackedArmsVec[i].lowerConfidenceBound
//                      << " ucb = " << hackedArmsVec[i].upperConfidenceBound
//                      << " total pulls="
//                      << hackedArmsVec[i].numberOfPulls << std::endl;
//        }
//
//        std::cout << "average pull " << UCB1.globalNumberOfPulls / armsVec.size() << std::endl;
//        std::cout << "sigma " << UCB1.globalSigma << std::endl;
//        std::cout << "best arm's estimate " << UCB1.arms.top().estimateOfMean << std::endl;
//        std::cout << UCB1.arms.top().id << std::endl;
//#endif

    }
    return 0;
}
