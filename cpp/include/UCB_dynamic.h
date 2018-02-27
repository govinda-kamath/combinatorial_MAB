//
// Created by Govinda Kamath on 2/26/18.
//

#ifndef COMBINATORIAL_MAB_UCB_DYNAMIC_H
#define COMBINATORIAL_MAB_UCB_DYNAMIC_H
#include <boost/heap/fibonacci_heap.hpp>

template <class templateArm>
struct compare_arms
{
    bool operator()(const templateArm &l, const templateArm &r) const
    {
        return r > l;
    }
};

template <class templateArm>
class UCBDynamic{
    /* UCB for the general case*/
public:
    unsigned numberOfArms;
    float logDeltaInverse;

    float globalSigma;
    unsigned long long  globalNumberOfPulls;
    float globalSumOfPulls;
    float globalSumOfSquaresOfPulls;
    unsigned numberOfBestArms;

    typedef typename boost::heap::fibonacci_heap<templateArm, boost::heap::compare< compare_arms<templateArm> > > fib_heap_ucb;

    std::vector<templateArm> armsContainer;
    fib_heap_ucb arms;
    std::vector<templateArm> topKArms;
    std::unordered_map < unsigned long, typename fib_heap_ucb::handle_type > handlesVec;





    UCBDynamic(std::vector<templateArm> &armsVec, float delta, unsigned nOfBestArms){
        armsContainer = armsVec;
        numberOfArms = armsContainer.size();
        logDeltaInverse = std::log(1/delta);
        globalSigma = NAN;
        globalNumberOfPulls = 0;
        globalSumOfPulls = 0;
        globalSumOfSquaresOfPulls = 0;
        numberOfBestArms = nOfBestArms;
    }

    void updateGlobalSigma(){
        globalSigma = std::sqrt((globalSumOfSquaresOfPulls/globalNumberOfPulls -
                                 std::pow(globalSumOfPulls/globalNumberOfPulls,2)));
    }

    void initialise(unsigned numberOfInitialPulls = 100){
        for (unsigned long armIndex = 0; armIndex< numberOfArms; armIndex++) {
#ifdef DEBUG_INIT
            if (armIndex%((int)(numberOfArms)/20) == 0){
                std::cout << "Initialized " << std::setprecision (15) << armIndex << " out of " << numberOfArms
                          << std::endl;
            }
#endif
            initialiseSingleArm( armsContainer[armIndex],  numberOfInitialPulls );
        }
        updateGlobalSigma();
#ifdef DEBUG_INIT
        std::cout << "Sigma after initialization " << globalSigma <<std::endl;
        std::cout << "Mean after initialization "
                  << globalSumOfPulls/globalNumberOfPulls <<std::endl;
        std::cout << "Mean Squared after initialization "
                  << std::pow(globalSumOfPulls/globalNumberOfPulls,2) <<std::endl;
        std::cout << "Second Moment after initialization "
                  << globalSumOfSquaresOfPulls/globalNumberOfPulls
                  << " No of pulls " << globalNumberOfPulls
                    <<std::endl;
#endif

        for (unsigned long armIndex = 0; armIndex < numberOfArms; armIndex++){
            addSingleArm(armsContainer[armIndex]);
        }
    }

    void initialiseSingleArm( templateArm singleArm, unsigned numberOfInitialPulls = 100){
        for (unsigned i = 0; i < numberOfInitialPulls; i++) {
            float observedSample(0);
            observedSample = singleArm.pullArm(0, globalNumberOfPulls, 0, false);
            globalSumOfPulls += observedSample;
            globalSumOfSquaresOfPulls += observedSample * observedSample;
        }
        globalNumberOfPulls += numberOfInitialPulls;

    }

    void addSingleArm( templateArm singleArm){
        singleArm.updateConfidenceIntervals(globalSigma, globalNumberOfPulls, logDeltaInverse);
        unsigned long armID;
        armID = singleArm.id;
        handlesVec[armID] = arms.push(singleArm);
    }


    void initialiseAndAddNewArm( templateArm newArm, unsigned numberOfInitialPulls = 100){
        initialiseSingleArm(newArm, numberOfInitialPulls);
        updateGlobalSigma();
        addSingleArm(newArm);

    }


    void runUCB(unsigned long maxIterations){
        unsigned bestArmCount = 0;
        unsigned long i(0);
        for (; i < maxIterations; i++){

//#define DEBUG_RUN
#ifdef DEBUG_RUN
            if (i%((int)maxIterations/200) == 0){
                templateArm bestArm = arms.top();
                arms.pop();
                templateArm secondBestArm = arms.top();
                arms.push(bestArm);

                float UCBofBestArm, LCBofBestArm;
                UCBofBestArm = bestArm.upperConfidenceBound;
                LCBofBestArm = bestArm.lowerConfidenceBound;

                std::cout << "NumberOfPulls " << globalNumberOfPulls << " out of " << maxIterations
                          << ". Best arm = " << bestArm.id
                          << ". Best arm UCB = " << UCBofBestArm
                          << ". LCB of second best arm  = " << LCBofBestArm
                          << ". No. pulls  = " <<  bestArm.numberOfPulls
                          << ". Estimate  = " <<  bestArm.estimateOfMean
                          << ". globalSigma = " << globalSigma
                          << std::endl;
            }
#endif
            bool bestArmFound = iterationOfUCB();
            if (bestArmFound){
                topKArms.push_back(arms.top());
                arms.pop();
                bestArmCount++;

#ifdef DEBUG_RUN
                std::cout << "Best arm number " << bestArmCount  << " Position " << i <<std::endl;
#endif
                if (bestArmCount==numberOfBestArms)
                    break;
            }
        }
        if (bestArmCount!=numberOfBestArms){
            std::cout<< "UCB Stopped before reaching optimal" << std::endl;
        }
#ifdef DEBUG_RUN
        std::cout << "Best arm number " << bestArmCount << " Position" << i
                  << " Max iter" << maxIterations  << std::endl;
#endif
        storeExtraTopArms(); //Storing extra arms
    }

    bool iterationOfUCB(){
        /*An iteration of UCB*/
        auto iter = arms.ordered_begin();
        unsigned long bestArmId = ((templateArm) (*iter)).id;//typecasting
        iter++;
        unsigned long secondBestArmId = ((templateArm) (*iter)).id;//typecasting

        auto handleBestArm = handlesVec[bestArmId];
        auto handleSecondBestArm = handlesVec[secondBestArmId];

//        auto handleBestArm = fib_heap_ucb::s_handle_from_iterator(iter);
//        iter++;
//        auto handleSecondBestArm = fib_heap_ucb::s_handle_from_iterator(iter);

        float UCBofBestArm, LCBofSecondBestArm;
        UCBofBestArm = (*handleBestArm).upperConfidenceBound;
        LCBofSecondBestArm = (*handleSecondBestArm).lowerConfidenceBound;

        if (UCBofBestArm < LCBofSecondBestArm){
            // Evaluating true mean of best arm
            std::unordered_map<std::string, float> result = (*handleBestArm).trueMeanUpdate();

            globalNumberOfPulls += result["effectiveDimension"];
            globalSumOfPulls += result["sumOfPulls"];
            globalSumOfSquaresOfPulls += result["sumOfSquaresPulls"];

            if ((*handleBestArm).estimateOfMean > LCBofSecondBestArm){
                std::cout<< "False trigger by " << (*handleBestArm).id << std::endl;
                arms.update(handleBestArm);
                return false;
            }

#ifdef DEBUG_RUN
            std::cout << "stopping UCB "<< std::setprecision (15)<< UCBofBestArm << "id " << bestArm.id <<  std::endl;
            std::cout << "stopping LCB " <<  LCBofSecondBestArm << "id " << secondBestArm.id << std::endl;
            std::cout << "best estimate"<< std::setprecision (15)<< bestArm.estimateOfMean << std::endl;
#endif
            return true;
        }else {
            float sample;
            sample = (*handleBestArm).pullArm(globalSigma, globalNumberOfPulls, logDeltaInverse, true);
            globalSumOfPulls += sample;
            globalSumOfSquaresOfPulls += std::pow(sample, 2);
            globalNumberOfPulls++;
            globalSigma = std::sqrt((globalSumOfSquaresOfPulls / globalNumberOfPulls -
                                     std::pow(globalSumOfPulls / globalNumberOfPulls, 2)));
            arms.update(handleBestArm);
            return  false;
        }
    }


    void storeExtraTopArms(){
        for (unsigned i=0; i < std::min(numberOfBestArms*5, numberOfArms - numberOfBestArms); i++) {
            topKArms.push_back(arms.top());
            arms.pop();
        }
    }

    // Never call this for large dataset!
    templateArm bestArm(){
        unsigned long bIndex = 0;
        float minTrueMean = armsContainer[0].trueMean();
        assert(("Dataset too large", armsContainer.size() < 1000));
        for (unsigned long i(0); i< armsContainer.size(); i++){
            float tmpTrueMean = armsContainer[i].trueMean();
            if ( tmpTrueMean < minTrueMean){
                minTrueMean = tmpTrueMean;
                bIndex = i;
            }
        }
        return armsContainer[bIndex];
    }
};

#endif //COMBINATORIAL_MAB_UCB_DYNAMIC_H
