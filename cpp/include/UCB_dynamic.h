//
// Created by Govinda Kamath on 2/28/18.
//

#ifndef COMBINATORIAL_MAB_UCB_DYNAMIC_H
#define COMBINATORIAL_MAB_UCB_DYNAMIC_H

#include <unordered_set>
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
    unsigned numberOfExtraArms;

    std::vector<templateArm> armsContainer;
    std::priority_queue<templateArm, std::vector<templateArm>, std::greater<templateArm> > arms;
    std::vector<templateArm> topKArms;
    std::unordered_set <unsigned long > armsToKeep;

    UCBDynamic(std::vector<templateArm> &armsVec, float delta, unsigned nOfBestArms, unsigned nOfExtraArms){
        armsContainer = armsVec;
        numberOfArms = armsContainer.size();
        logDeltaInverse = std::log(1/delta);
        globalSigma = NAN;
        globalNumberOfPulls = 0;
        globalSumOfPulls = 0;
        globalSumOfSquaresOfPulls = 0;
        numberOfBestArms = nOfBestArms;
        numberOfExtraArms = nOfExtraArms;


    }

    void updateGlobalSigma(){
        globalSigma = std::sqrt((globalSumOfSquaresOfPulls/globalNumberOfPulls -
                                 std::pow(globalSumOfPulls/globalNumberOfPulls,2)));
    }

    void initialise(unsigned numberOfInitialPulls = 100){
        for (unsigned long armIndex = 0; armIndex< numberOfArms; armIndex++) {
            initialiseSingleArm( armsContainer[armIndex],  numberOfInitialPulls );

//#define DEBUG_INIT
#ifdef DEBUG_INIT
            if (armIndex%((int)(numberOfArms)/20) == 0){
                std::cout << "Initialized " << std::setprecision (15) << armIndex << " out of " << numberOfArms
                          << " Value: " << armsContainer[armIndex].estimateOfMean << std::endl;
            }

#endif
        }
        updateGlobalSigma();

        for (unsigned long armIndex = 0; armIndex < numberOfArms; armIndex++){
            addSingleArm(armsContainer[armIndex]);
        }
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
    }

    void initialiseSingleArm( templateArm &singleArm, unsigned numberOfInitialPulls = 100){
        for (unsigned i = 0; i < numberOfInitialPulls; i++) {
            float observedSample(0);
            observedSample = singleArm.pullArm(0, globalNumberOfPulls, 0, false);
            globalSumOfPulls += observedSample;
            globalSumOfSquaresOfPulls += observedSample * observedSample;
        }
        globalNumberOfPulls += numberOfInitialPulls;

    }

    void addSingleArm( templateArm &singleArm){
        singleArm.updateConfidenceIntervals(globalSigma, globalNumberOfPulls, logDeltaInverse);
        armsToKeep.insert(singleArm.id);
        arms.push(singleArm);
    }


    void initialiseAndAddNewArm( templateArm &newArm, unsigned numberOfInitialPulls = 100){
        initialiseSingleArm(newArm, numberOfInitialPulls);
        updateGlobalSigma();
        addSingleArm(newArm);

    }

    void markForRemoval(unsigned long armID){
        auto search = armsToKeep.find(armID);
        if (search != armsToKeep.end())
            armsToKeep.erase(armID);
    }

    templateArm topValidArm(){
        bool topValidArmFound(false);
        do {
            unsigned  long topArmID = arms.top().id;
            if (armsToKeep.find(topArmID) != armsToKeep.end())
                topValidArmFound=true;
            else{
                arms.pop();
            }
        } while((!topValidArmFound) or (arms.empty()) );
        if (arms.empty()){
            throw std::runtime_error("[Unexpected behaviour]: Arms Priority Queue empty.");
        }
        return arms.top();
    }

    void runUCB(unsigned long maxIterations){
        unsigned bestArmCount = 0;
        unsigned long i(0);
        for (; i < maxIterations; i++){

//#define DEBUG_RUN
#ifdef DEBUG_RUN
            if (i%((int)maxIterations/200) == 0){
                templateArm bestArm = topValidArm();
                arms.pop();
                templateArm secondBestArm = topValidArm();
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
                templateArm topArm = topValidArm();
                topKArms.push_back(topArm);
                arms.pop();
                markForRemoval(topArm.id);
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
        templateArm bestArm = topValidArm();
        arms.pop();
        templateArm secondBestArm = topValidArm();
        float UCBofBestArm, LCBofSecondBestArm;
        UCBofBestArm = bestArm.upperConfidenceBound;
        LCBofSecondBestArm = secondBestArm.lowerConfidenceBound;

        if (UCBofBestArm < LCBofSecondBestArm){
            // Evaluating true mean of best arm
            std::unordered_map<std::string, float> result = bestArm.trueMeanUpdate();
            bestArm.estimateOfMean = result["sumOfPulls"]/result["effectiveDimension"];
            bestArm.upperConfidenceBound = bestArm.estimateOfMean;
            bestArm.lowerConfidenceBound = bestArm.estimateOfMean;
            globalNumberOfPulls += result["effectiveDimension"];
            globalSumOfPulls += result["sumOfPulls"];
            globalSumOfSquaresOfPulls += result["sumOfSquaresPulls"];

            if (bestArm.estimateOfMean > LCBofSecondBestArm){
                std::cout<< "False trigger by " << bestArm.id << std::endl;
                arms.push(bestArm);
                return false;
            }

            arms.push(bestArm);
#ifdef DEBUG_RUN
            std::cout << "stopping UCB "<< std::setprecision (15)<< UCBofBestArm << "id " << bestArm.id <<  std::endl;
            std::cout << "stopping LCB " <<  LCBofSecondBestArm << "id " << secondBestArm.id << std::endl;
            std::cout << "best estimate"<< std::setprecision (15)<< bestArm.estimateOfMean << std::endl;
#endif
            return true;
        }else {
            float sample;
            sample = bestArm.pullArm(globalSigma, globalNumberOfPulls, logDeltaInverse, true);
            globalSumOfPulls += sample;
            globalSumOfSquaresOfPulls += std::pow(sample, 2);
            globalNumberOfPulls++;
            globalSigma = std::sqrt((globalSumOfSquaresOfPulls / globalNumberOfPulls -
                                     std::pow(globalSumOfPulls / globalNumberOfPulls, 2)));
            arms.push(bestArm);
            return  false;
        }
    }


    void storeExtraTopArms(){
        for (unsigned i=0; i < numberOfExtraArms; i++) {
            topKArms.push_back(topValidArm());
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
