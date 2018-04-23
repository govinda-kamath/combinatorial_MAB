//
// Created by Govinda Kamath on 2/28/18.
//

#ifndef COMBINATORIAL_MAB_UCB_DYNAMIC_H
#define COMBINATORIAL_MAB_UCB_DYNAMIC_H

#include "utils.h"
#include <unordered_set>
#include <unordered_map>
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
    unsigned sampleSize;

    std::vector<templateArm> armsContainer;
    std::priority_queue<templateArm, std::vector<templateArm>, std::greater<templateArm> > arms;
    std::priority_queue<templateArm, std::vector<templateArm>, std::greater<templateArm> > armsBrute;
    std::vector<templateArm> topKArms;
    std::unordered_set <unsigned long > armsToKeep;
    std::unordered_map<unsigned long, utils::ArmConditions> armStates;


    // Main Constructor
    UCBDynamic(std::vector<templateArm> &armsVec, float delta, unsigned nOfBestArms, unsigned nOfExtraArms, unsigned sSize){
        armsContainer = armsVec;
        numberOfArms = armsContainer.size();
        logDeltaInverse = std::log(1/delta);
        globalSigma = NAN;
        globalNumberOfPulls = 0;
        globalSumOfPulls = 0;
        globalSumOfSquaresOfPulls = 0;
        numberOfBestArms = nOfBestArms;
        numberOfExtraArms = nOfExtraArms;
        sampleSize = sSize;
        armStates.reserve(armsContainer.size()*2);
    }

    void updateGlobalSigma(){
        globalSigma = std::sqrt((globalSumOfSquaresOfPulls/globalNumberOfPulls -
                                 std::pow(globalSumOfPulls/globalNumberOfPulls,2)));
    }

    // Step 1 of UCB
    void initialise(unsigned numberOfInitialPulls = 100){
        for (unsigned long armIndex = 0; armIndex< numberOfArms; armIndex++) {
            if (armIndex%10000==0)
                std::cout << armIndex << "\t";
            initialiseSingleArm( armsContainer[armIndex],  numberOfInitialPulls );
        }
        updateGlobalSigma();
        std::cout << "Global sigma after initialization =  " << globalSigma << std::endl;
        std::cout << "Global Number Of Pulls =  " << globalNumberOfPulls << std::endl;
        for (unsigned long armIndex = 0; armIndex < numberOfArms; armIndex++){
            addSingleArm(armsContainer[armIndex]);
        }
    }

    // Used by Step 1 of UCB
    void initialiseSingleArm( templateArm &singleArm, unsigned numberOfInitialPulls = 100){
        std::pair<float, float> sample;
        sample = singleArm.pullArm(0, NAN, 0, false, numberOfInitialPulls);
        globalSumOfPulls += sample.first;
        globalSumOfSquaresOfPulls += sample.second;
        globalNumberOfPulls += numberOfInitialPulls;

    }

    // Used by Step 1 of UCB
    void addSingleArm( templateArm &singleArm){
        singleArm.updateConfidenceIntervals(globalSigma, globalNumberOfPulls, logDeltaInverse);
        armStates[singleArm.id] = utils::ArmConditions(singleArm.numberOfPulls, singleArm.sumOfPulls,
                                         singleArm.sumOfSquaresOfPulls, singleArm.trueMeanValue);
        armsToKeep.insert(singleArm.id);
        arms.push(singleArm);
    }

    // Dynamic part of UCB. Uaed to add new arm
    void initialiseAndAddNewArm( templateArm &newArm, unsigned numberOfInitialPulls = 100){
        initialiseSingleArm(newArm, numberOfInitialPulls);
        newArm.updateConfidenceIntervals(globalSigma, globalNumberOfPulls, logDeltaInverse);
        armStates[newArm.id] = utils::ArmConditions(newArm.numberOfPulls, newArm.sumOfPulls,
                                           newArm.sumOfSquaresOfPulls, newArm.trueMeanValue);
        armsToKeep.insert(newArm.id);
        arms.push(newArm);
    }

    // Dynamic part of UCB. Removes existing arms
    void markForRemoval(unsigned long armID){
        auto search = armsToKeep.find(armID);
        if (search != armsToKeep.end())
            armsToKeep.erase(armID);
    }

    // Dynamic part of UCB. Returns top arm (which has not been removed)
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

    // Step 2 of UCB : Main step
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
            // Run a "single" step of UCB
            bool bestArmFound = iterationOfUCB();

            // Update if best arm is found
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

        // Debug only
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

        if (UCBofBestArm == NAN){
            std::cout << "Damn Nan" << bestArm.id << std::endl;
        }
        if (UCBofBestArm < LCBofSecondBestArm){

            /* Evaluating true mean of best arm
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
            */
#ifdef DEBUG_RUN
            std::cout << "stopping UCB "<< std::setprecision (15)<< UCBofBestArm << "id " << bestArm.id <<  std::endl;
            std::cout << "stopping LCB " <<  LCBofSecondBestArm << "id " << secondBestArm.id << std::endl;
            std::cout << "best estimate"<< std::setprecision (15)<< bestArm.estimateOfMean << std::endl;
#endif

            arms.push(bestArm);
            return true;
        }
        // Run a iteration if the top arm is not the best arm
        else {
            std::pair<float, float> sample;
            sample = bestArm.pullArm(globalSigma, globalNumberOfPulls, logDeltaInverse, true, sampleSize, LCBofSecondBestArm);

            globalSumOfPulls += sample.first;
            globalSumOfSquaresOfPulls += sample.second;
            globalNumberOfPulls += sampleSize;
            float a, b, c;
            a = globalSumOfSquaresOfPulls / globalNumberOfPulls;
            b = globalSumOfPulls / globalNumberOfPulls;
            c = std::pow(b, 2);
            globalSigma = std::sqrt(( a - c ));
            if (globalSigma == NAN){
                std::cout << "second point is NAN" << std::endl;
            }
            arms.push(bestArm);
            //Update arm status
            unsigned long numArmPulls = bestArm.numberOfPulls;
            float armSumOfPulls = bestArm.sumOfPulls;
            float armSumOfSquaresOfPulls = bestArm.sumOfSquaresOfPulls;
            if(armStates.find(bestArm.id)!= armStates.end()){
                armStates[bestArm.id].update(numArmPulls, armSumOfPulls,armSumOfSquaresOfPulls);
            }
            else{
                throw std::runtime_error("[Unexpected behaviour]: Best arm's state not found.");
            }

            return  false;
        }
    }


    // Stores extra top arms for evaluations. These arms are not optimally chosen
    void storeExtraTopArms(){
        for (unsigned i=0; i < numberOfExtraArms; i++) {
            topKArms.push_back(topValidArm());
            arms.pop();
        }
    }


    /*
     *  Debugging: Brute force methods.
     */
    // For brute force only
    void initialiseAndAddNewArmBrute( templateArm &newArm, unsigned numberOfInitialPulls = 100){
        initialiseAndAddNewArm( newArm, numberOfInitialPulls);
        float tmpTrueMean = newArm.trueMean();
        armsBrute.push(newArm);
    }

    // For brute force only Initial Step
    void armsKeepFromArmsContainerBrute(){
        for (unsigned long index(0); index < armsContainer.size(); index++){
            float tmpTrueMean = armsContainer[index].trueMean();
            armsContainer[index].lowerConfidenceBound = tmpTrueMean;
            armsContainer[index].upperConfidenceBound = tmpTrueMean;
            armsContainer[index].estimateOfMean = tmpTrueMean;
            armsContainer[index].estimateOfSecondMoment = tmpTrueMean*tmpTrueMean;
            armsToKeep.insert(armsContainer[index].id);
            armStates[armsContainer[index].id] = utils::ArmConditions(1, tmpTrueMean, tmpTrueMean*tmpTrueMean,
                    armsContainer[index].trueMeanValue);
            armsBrute.push(armsContainer[index]);

        }
    }

    // For brute force only
    templateArm topValidArmBrute(){
        bool topValidArmFound(false);
        do {
            unsigned  long topArmID = armsBrute.top().id;
            if (armsToKeep.find(topArmID) != armsToKeep.end())
                topValidArmFound=true;
            else{
                armsBrute.pop();
            }
        } while((!topValidArmFound) or (armsBrute.empty()) );
        if (armsBrute.empty()){
            throw std::runtime_error("[Unexpected behaviour]: Arms Priority Queue empty.");
        }
        return armsBrute.top();
    }

    // For brute force only
    templateArm bruteBestArm(){
        templateArm bestArm = topValidArmBrute();
        armsBrute.pop();
        return bestArm;
    }
};

#endif //COMBINATORIAL_MAB_UCB_DYNAMIC_H
