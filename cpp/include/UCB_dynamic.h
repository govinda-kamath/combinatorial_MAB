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
    unsigned long  globalNumberOfPulls;
    float globalSumOfPulls;
    float globalSumOfSquaresOfPulls;
    unsigned numberOfBestArms;
    unsigned numberOfExtraArms;
    unsigned sampleSize;

    std::vector<templateArm> armsContainer;
    std::priority_queue<templateArm, std::vector<templateArm>, std::greater<templateArm> > arms;
    std::unordered_map<unsigned long, unsigned long> finalNumberOfPulls;
    std::vector<unsigned long> finalSortedOrder;

    std::priority_queue<templateArm, std::vector<templateArm>, std::greater<templateArm> > armsBrute;
    std::vector<templateArm> topKArms;
    std::unordered_set <unsigned long > armsToKeep;
    std::unordered_map<unsigned long, utils::ArmConditions> armStates;
    std::unordered_map<unsigned long, long int> finalNumberOfArmPulls;


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
    void initialise (unsigned numberOfInitialPulls = 100){
        for (unsigned long armIndex = 0; armIndex< numberOfArms; armIndex++) {
            initialiseSingleArm( armsContainer[armIndex],  numberOfInitialPulls );
        }
        updateGlobalSigma();
//        std::cout << "Global sigma after initialization =  " << globalSigma << std::endl;
//        std::cout << "Global Number Of Pulls =  " << globalNumberOfPulls << std::endl;
//        std::cout << "Adding to container ";
        for (unsigned long armIndex = 0; armIndex < numberOfArms; armIndex++){
//            std::cout<< armIndex << " ";
            addSingleArm(armsContainer[armIndex]);
        }
//        std::cout<<std::endl;
    }

    // Used by Step 1 of UCB
    void initialiseSingleArm( templateArm &singleArm, unsigned numberOfInitialPulls = 100){
        std::pair<float, float> sample;
        sample = singleArm.pullArm(0, NAN, logDeltaInverse, false, numberOfInitialPulls, -1); // Huge change made here
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

    // Dynamic part of UCB. Used to add new arm
    void initialiseAndAddNewArm( templateArm &newArm, unsigned numberOfInitialPulls = 100){
        if (numberOfInitialPulls>0){
            initialiseSingleArm(newArm, numberOfInitialPulls);
            newArm.updateConfidenceIntervals(globalSigma, globalNumberOfPulls, logDeltaInverse);
        }
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
        else{
            std::cout << "Not able to erase armID " << armID << std::endl;
        }
    }

    // Dynamic part of UCB. Returns top arm (which has not been removed)
    templateArm topValidArm(){
        bool topValidArmFound(false);
        do {
            unsigned  long topArmID = arms.top().id;
            if (armsToKeep.find(topArmID) != armsToKeep.end())
                topValidArmFound=true;
            else{
                finalNumberOfPulls[arms.top().id] = arms.top().numberOfPulls;
                arms.pop();
            }
        } while((!topValidArmFound) and (!arms.empty()) );
        if (arms.empty()){
            std::cout << "Arms Priority Queue empty." << std::endl;
            templateArm tmpArm(-1);
            return tmpArm;
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
//                std::cout << "best arm " << topArm.id << std::endl;

                // Storing the results
                topKArms.push_back(topArm);
                finalSortedOrder.push_back(topArm.id);
                finalNumberOfPulls[topArm.id] = topArm.numberOfPulls;

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
            std::cout<< "UCB Stopped before reaching optimal. Adding the top arms anyways" << std::endl;
            for (unsigned int remain = 0; remain<numberOfBestArms-bestArmCount; remain++){
                templateArm topArm = topValidArm();
                topKArms.push_back(topArm);
                finalSortedOrder.push_back(topArm.id);
                finalNumberOfPulls[topArm.id] = topArm.numberOfPulls;
                arms.pop();
                markForRemoval(topArm.id);
            }
        }
        std::cout << "Ran for " << i << " steps" << std::endl;
#ifdef DEBUG_RUN
        std::cout << "Best arm number " << bestArmCount << " Position" << i
                  << " Max iter" << maxIterations  << std::endl;
#endif
    }


    bool iterationOfUCB(){
        /*An iteration of UCB*/
        templateArm bestArm = topValidArm();
        arms.pop();
        templateArm secondBestArm = topValidArm();
        float UCBofBestArm, LCBofSecondBestArm;
        UCBofBestArm = bestArm.upperConfidenceBound;
        LCBofSecondBestArm = secondBestArm.lowerConfidenceBound;
        float del = UCBofBestArm-LCBofSecondBestArm;
//        if ( globalNumberOfPulls%23==1 ) {
//            std::cout << "POW! \n\t"
//                      << bestArm.id << " "
//                      << bestArm.lowerConfidenceBound << " "
//                      << bestArm.estimateOfMean << " "
//                      << bestArm.upperConfidenceBound << " "
//                      << bestArm.numberOfPulls << "\n\t"
//                      << secondBestArm.id << " "
//                      << secondBestArm.lowerConfidenceBound << " "
//                      << secondBestArm.estimateOfMean << " "
//                      << secondBestArm.upperConfidenceBound << " "
//                      << secondBestArm.numberOfPulls << " "
//                      << globalNumberOfPulls / arms.size()
//                      << " " << del << std::endl;
//        }
        if (del == NAN){
            std::cout << "Damn the NAN" << bestArm.id << std::endl;
        }
        if (UCBofBestArm < LCBofSecondBestArm){
//            std::cout << "POW! \n\t"
//                      << bestArm.id << " "
//                      << bestArm.lowerConfidenceBound << " "
//                      << bestArm.estimateOfMean << " "
//                      << bestArm.upperConfidenceBound << " "
//                      << bestArm.numberOfPulls << "\n\t"
//                      << secondBestArm.id << " "
//                      << secondBestArm.lowerConfidenceBound << " "
//                      << secondBestArm.estimateOfMean << " "
//                      << secondBestArm.upperConfidenceBound << " "
//                      << secondBestArm.numberOfPulls << " "
//                      << globalNumberOfPulls / arms.size()
//                      << " " << del << std::endl;
//            std::cout << "BEST ARM ARM!!"
//                    ""
//                      << "\n" << bestArm.id << " "
//                      << secondBestArm.id << " "
//                      << UCBofBestArm << " "
//                      << LCBofSecondBestArm << "\t"
//                      << del <<std::endl;
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
            if(sample.first!=0) {
                globalSumOfPulls += sample.first;
                globalSumOfSquaresOfPulls += sample.second;
                globalNumberOfPulls += sampleSize;
            }
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
        for (unsigned i=0; arms.size()>0; i++) {
            templateArm topArm = topValidArm();
            if (topArm.id == -1){
                std::cout << "Breaking Bad Here" << std::endl;
                break;
            }
            finalSortedOrder.push_back(topArm.id);
            finalNumberOfPulls[topArm.id] = topArm.numberOfPulls;
            if (i < numberOfExtraArms)
                topKArms.push_back(topArm);
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
            if (armsToKeep.find(topArmID) != armsToKeep.end()){
                topValidArmFound=true;
            }
            else{
                armsBrute.pop();
            }
        } while((!topValidArmFound) and (!armsBrute.empty()) );
        if (armsBrute.empty()){
            throw std::runtime_error("[Unexpected behaviour]: Arms Priority Queue empty.");
        }
        return armsBrute.top();
    }

    // For brute force only
    std::vector<templateArm> bruteBestArms(){
        std::vector<templateArm> topArms;
        for(unsigned long i = 0 ; i < numberOfBestArms+numberOfExtraArms; i++){
            templateArm topArm = topValidArmBrute();
            armsBrute.pop();
            topArms.push_back(topArm);
            markForRemoval(topArm.id);

        }
//        std::cout<<"Returning " << topArms.size()<<  std::endl;
        return topArms;
    }
};

#endif //COMBINATORIAL_MAB_UCB_DYNAMIC_H
