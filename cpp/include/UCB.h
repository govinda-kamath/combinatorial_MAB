//
// Created by Vivek Kumar Bagaria on 1/31/18.
//

#ifndef COMBINATORIAL_MAB_UCB_H
#define COMBINATORIAL_MAB_UCB_H


template <class templateArm>
class UCB{
    /* UCB for the general case*/
public:
    unsigned numberOfArms;
    float logDeltaInverse;

    float globalSigma;
    unsigned long long  globalNumberOfPulls;
    float globalSumOfPulls;
    float globalSumOfSquaresOfPulls;
    unsigned numberOfBestArms;

    std::vector<templateArm> armsContainer;
    std::priority_queue<templateArm, std::vector<templateArm>, std::greater<templateArm> > arms;
    std::vector<templateArm> topKArms;
    std::mutex initializeMutex;
    UCB(std::vector<templateArm> &armsVec, float delta, unsigned nOfBestArms){

        armsContainer = armsVec;
        numberOfArms = armsContainer.size();

        logDeltaInverse = std::log(1/delta);

        globalSigma = NAN;
        globalNumberOfPulls = 0;
        globalSumOfPulls = 0;
        globalSumOfSquaresOfPulls = 0;
        numberOfBestArms = nOfBestArms;
    }

    void initialiseFewArm(unsigned long armIndexStart, unsigned long armIndexEnd, unsigned numberOfInitialPulls){

        float localSumOfPulls, localSumOfSquaresOfPulls;
        localSumOfPulls = 0;
        localSumOfSquaresOfPulls = 0;
        // Pulling an arm numberOfInitialPulls times
        for (unsigned long armIndex = armIndexStart; armIndex< armIndexEnd; armIndex++) {
#ifdef DEBUG_INIT
            if (armIndex%((int)(armIndexEnd-armIndexStart)/20) == 0){
                std::cout << "Initialized " << std::setprecision (15) << armIndex << " out of " << armIndexEnd - armIndexStart
                          << std::endl;
            }
#endif
            for (unsigned i = 0; i < numberOfInitialPulls; i++) {
                float observedSample(0);
                observedSample = armsContainer[armIndex].pullArm(0, globalNumberOfPulls, 0, false);
                localSumOfPulls += observedSample;
                localSumOfSquaresOfPulls += observedSample * observedSample;

            }
        }
        // locking the global variables which have to be updated
        {
//            std::lock_guard<std::mutex> guard(initializeMutex);
            globalSumOfPulls += localSumOfPulls;
            globalSumOfSquaresOfPulls += localSumOfSquaresOfPulls;
        }
    }
    void initialise(unsigned numberOfInitialPulls = 100){

        unsigned numberOfThreads = 1 ;//std::thread::hardware_concurrency();
        std::vector<std::thread> initThreads(numberOfThreads);

        unsigned long chunkSize = (numberOfArms/numberOfThreads);
        for(unsigned t = 0; t < numberOfThreads; t++){
            unsigned long armIndexStart = t*chunkSize;
            unsigned long armIndexEnd = (t+1)*chunkSize;
            initThreads[t] = std::thread(&UCB::initialiseFewArm, this,
                                         armIndexStart, armIndexEnd, numberOfInitialPulls);
        }

        for(unsigned t = 0; t < numberOfThreads; t++){
            initThreads[t].join();
        }

        globalNumberOfPulls = numberOfInitialPulls*numberOfArms;
        globalSigma = std::sqrt((globalSumOfSquaresOfPulls/globalNumberOfPulls -
                                 std::pow(globalSumOfPulls/globalNumberOfPulls,2)));
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

        for (unsigned long index = 0; index < numberOfArms; index++){
            armsContainer[index].updateConfidenceIntervals(globalSigma, globalNumberOfPulls, logDeltaInverse);
            arms.push(armsContainer[index]);
        }

    }


    void runUCB(unsigned long maxIterations){
        unsigned bestArmCount = 0;
        unsigned long i(0);
        for (; i < maxIterations; i++){

            if (i%((int)maxIterations/200) == 0){
                templateArm bestArm = arms.top();
                arms.pop();
                templateArm secondBestArm = arms.top();
                arms.push(bestArm);

                float UCBofBestArm, LCBofBestArm;
                UCBofBestArm = bestArm.upperConfidenceBound;
                LCBofBestArm = bestArm.lowerConfidenceBound;
//#define DEBUG_RUN
#ifdef DEBUG_RUN
                std::cout << "NumberOfPulls " << globalNumberOfPulls << " out of " << maxIterations
                          << ". Best arm = " << bestArm.id
                          << ". Best arm UCB = " << UCBofBestArm
                          << ". LCB of second best arm  = " << LCBofBestArm
                          << ". No. pulls  = " <<  bestArm.numberOfPulls
                          << ". Estimate  = " <<  bestArm.estimateOfMean
                          << ". globalSigma = " << globalSigma
                          << std::endl;
#endif
            }

            bool bestArmFound;
            bestArmFound = iterationOfUCB();
            if (bestArmFound){
                topKArms.push_back(arms.top());
                arms.pop();
                bestArmCount++;

#ifdef DEBUG_RUN
                std::cout << "Best arm number " << bestArmCount
                          << " Position " << i <<std::endl;
#endif
                if (bestArmCount==numberOfBestArms)
                    break;
            }
        }
        if (bestArmCount!=numberOfBestArms){
            std::cout<< "UCB Stopped before reaching optimal" << std::endl;
        }
#ifdef DEBUG_RUN
        std::cout << "Best arm number "
                  << bestArmCount << " Position" << i
                  << " Max iter" << maxIterations
                  << std::endl;
#endif
        storeExtraTopArms(); //Storing extra arms
    }

    void storeExtraTopArms(){
        for (unsigned i=0; i < std::min(numberOfBestArms*5, numberOfArms - numberOfBestArms); i++) {
            topKArms.push_back(arms.top());
            arms.pop();
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

#ifdef DEBUG_RUN
            std::cout << "stopping UCB "<< std::setprecision (15)<< UCBofBestArm << "id " << bestArm.id <<  std::endl;
            std::cout << "stopping LCB " <<  LCBofSecondBestArm << "id " << secondBestArm.id << std::endl;
#endif

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
            std::cout << "best estimate"<< std::setprecision (15)<< bestArm.estimateOfMean << std::endl;
#endif
            return true;
        }else {
            float sample;
            sample = bestArm.pullArm(globalSigma, globalNumberOfPulls, logDeltaInverse, true);
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

    templateArm bestArm(){
        unsigned long bIndex = 0;
        float minTrueMean = armsContainer[0].trueMean();
        for (unsigned long i(0); i< armsContainer.size(); i++){
            float tmpTrueMean = armsContainer[i].trueMean();
//            std::cout << i << " " << tmpTrueMean << std::endl;
            if ( tmpTrueMean < minTrueMean){
                minTrueMean = tmpTrueMean;
                bIndex = i;
            }
        }
        return armsContainer[bIndex];
    }
};


#endif //COMBINATORIAL_MAB_UCB_H
