#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cassert>
#include <queue>

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
    Point(std::vector <float> p){
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
    SquaredEuclideanPoint(std::vector<float> p) : Point(p){}

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
    }

    ~Arm(){
//        delete &point; //no clue why this errors out
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

        sample = point->sampledDistance(p1);

        numberOfPulls ++;
        sumOfPulls += sample;
        estimateOfMean = sumOfPulls/numberOfPulls;
        SumOfSquaresOfPulls += sample*sample;
        estimateOfSecondMoment = SumOfSquaresOfPulls/numberOfPulls;

        if (update)
            updateConfidenceIntervals(globalSigma,logDeltaInverse);

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

    UCB(std::vector<templateArm> &armsVec, float delta){

        armsContainer = armsVec;
        numberOfArms = armsContainer.size();

        logDeltaInverse = std::log(1/delta);

        globalSigma = NAN;
        globalNumberOfPulls = 0;
        globalSumOfPulls = 0;
        globalSumOfSquaresOfPulls = 0;

    }

    void initialise(int numberOfInitialPulls = 100){

        for (unsigned long index = 0; index < armsContainer.size(); index++){
            for (unsigned i = 0; i < numberOfInitialPulls; i++) {
                float observedSample(0);
                observedSample = armsContainer[index].pullArm(0, 0, false);
                globalSumOfPulls += observedSample;
                globalSumOfSquaresOfPulls += observedSample * observedSample;
            }
        }
        globalNumberOfPulls = numberOfInitialPulls*armsContainer.size();
        globalSigma = std::sqrt((globalSumOfSquaresOfPulls/globalNumberOfPulls -
                                    std::pow(globalSumOfPulls/globalNumberOfPulls,2)));
#ifdef DEBUG
        std::cout << "Sigma after initialization " << globalSigma <<std::endl;
        std::cout << "Mean after initialization " << globalSumOfPulls/globalNumberOfPulls <<std::endl;
        std::cout << "Mean Squared after initialization " << std::pow(globalSumOfPulls/globalNumberOfPulls,2) <<std::endl;
        std::cout << "Second Moment after initialization " << globalSumOfSquaresOfPulls/globalNumberOfPulls <<std::endl;
#endif

        for (unsigned long index = 0; index < armsContainer.size(); index++){

            armsContainer[index].updateConfidenceIntervals(globalSigma, logDeltaInverse);
            arms.push(armsContainer[index]);
        }

    }

    void runUCB(unsigned maxIterations){
        for (unsigned i(0); i < maxIterations; i++){
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
        float UCBofBestArm, LCBofSecondBestArm;
        UCBofBestArm = bestArm.upperConfidenceBound;
        LCBofSecondBestArm = secondBestArm.lowerConfidenceBound;
        if (UCBofBestArm < LCBofSecondBestArm){
            arms.push(bestArm);
            return true;
        } else {
            float sample;
            sample = bestArm.pullArm(globalSigma, logDeltaInverse);
            globalSumOfPulls += sample;
            globalSumOfSquaresOfPulls += std::pow(sample,2);
            globalNumberOfPulls++;
            globalSigma = std::sqrt((globalSumOfSquaresOfPulls/globalNumberOfPulls -
                                     std::pow(globalSumOfPulls/globalNumberOfPulls,2)));
            arms.push(bestArm);
        }
        return  false;
    }
};



int main(int argc, char *argv[]){

    std::cout << "We have entered " << argc
         << " arguments." << std::endl;



    std::string filePath(argv[1]), line;
    int numberOfInitialPulls(atoi(argv[2]));
    float delta(atof(argv[3]));

//    filePath = "/Users/govinda/Code/combinatorial_MAB/test_dataset/1000_images.txt";
//    filePath ="/data/MAB/work/dataset/test_dataset/basic_io_dataset/10k_images.txt";
    std::fstream fileReader(filePath.c_str());
    unsigned long pointIndex(0);

    std::vector<SquaredEuclideanPoint > pointsVec;
    std::vector<ArmKNN<SquaredEuclideanPoint> > armsVec;


    while(std::getline(fileReader, line)){
        float tmpValue;

        std::vector<float> tmpVec;
        std::stringstream ss(line);
        while (ss >> tmpValue){
            tmpVec.push_back(tmpValue);
        }
        SquaredEuclideanPoint tmpPoint(tmpVec);

        pointsVec.push_back(tmpPoint);
        pointIndex++;

    }

    for (unsigned i(1); i < pointsVec.size(); i++){
        ArmKNN<SquaredEuclideanPoint> tmpArm(i-1, pointsVec[i], pointsVec[0]);
        armsVec.push_back(tmpArm);
    }


    UCB<ArmKNN<SquaredEuclideanPoint> > UCB1(armsVec,delta);

    UCB1.initialise(numberOfInitialPulls);
    std::vector<ArmKNN<SquaredEuclideanPoint> > &hackedArmsVec = Container(UCB1.arms);
#ifdef DEBUG
    for(unsigned i=0; i< hackedArmsVec.size(); i++){
        std::cout << hackedArmsVec[i].id <<" True Mean= "<< armsVec[hackedArmsVec[i].id].trueMean()
                  << " sigma = " << std::sqrt((hackedArmsVec[i].SumOfSquaresOfPulls/hackedArmsVec[i].numberOfPulls -
                                               std::pow(hackedArmsVec[i].sumOfPulls/hackedArmsVec[i].numberOfPulls,2)))
                  << " estimate = " << hackedArmsVec[i].estimateOfMean << " total pulls="
                  << hackedArmsVec[i].numberOfPulls << std::endl;
    }

    std::cout << "average pull " << UCB1.globalNumberOfPulls/armsVec.size()<<std::endl;
    std::cout << "sigma " << UCB1.globalSigma<<std::endl;
    std::cout << "best arm's estimate "<<UCB1.arms.top().estimateOfMean << std::endl;
    std::cout << UCB1.arms.top().id << std::endl;
#endif
    UCB1.runUCB(1000);

//    std::vector<ArmKNN<SquaredEuclideanPoint> > &hackedArmsVec = Container(UCB1.arms);
#ifdef DEBUG
    for(unsigned i=0; i< hackedArmsVec.size(); i++){
        std::cout << hackedArmsVec[i].id <<" True Mean= "<< armsVec[hackedArmsVec[i].id].trueMean()
                << " sigma = " << std::sqrt((hackedArmsVec[i].SumOfSquaresOfPulls/hackedArmsVec[i].numberOfPulls -
                                             std::pow(hackedArmsVec[i].sumOfPulls/hackedArmsVec[i].numberOfPulls,2)))
                << " estimate = " << hackedArmsVec[i].estimateOfMean << " total pulls="
                  << hackedArmsVec[i].numberOfPulls << std::endl;
    }

    std::cout << "average pull " << UCB1.globalNumberOfPulls/armsVec.size()<<std::endl;
    std::cout << "sigma " << UCB1.globalSigma<<std::endl;
    std::cout << "best arm's estimate "<<UCB1.arms.top().estimateOfMean << std::endl;
    std::cout << UCB1.arms.top().id << std::endl;
#endif
    return 0;
}
