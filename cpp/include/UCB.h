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
    float globalNumberOfPulls;
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

        std::cout << "num initial pulls "<< numberOfInitialPulls << std::endl;
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
    void initialise(unsigned numberOfInitialPulls = 100){

        unsigned numberOfThreads = 1 ;//std::thread::hardware_concurrency();
//        std::cout << numberOfThreads << " number of threads used in each batch of initialization.\n";
        std::vector<std::thread> initThreads(numberOfThreads);

        unsigned long chunkSize = (numberOfArms/numberOfThreads);
        for(unsigned t = 0; t < numberOfThreads; t++){
            unsigned long armIndexStart = t*chunkSize;
            unsigned long armIndexEnd = (t+1)*chunkSize;
            initThreads[t] = std::thread(&UCB::initialiseFewArm, this,
                                         armIndexStart, armIndexEnd, numberOfInitialPulls);
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
        unsigned bestArmCount = 0;
        unsigned long i(0);
        for (; i < maxIterations; i++){

//            if ( armsContainer[7374].lowerConfidenceBound < arms.top().lowerConfidenceBound ){
//                std::vector<templateArm> tmp;
//                std::cout << std::setprecision (15) << "i=" << i << std::endl;
//                std::cout << " p size" << arms.size() << std::endl;
//
//
//                for (unsigned j(0); j < 20000; j++){
//
//                    std::cout << " j = " << j
//                            << "UCB= " <<    arms.top().upperConfidenceBound
//                            << "LCB= " <<    arms.top().lowerConfidenceBound
//                              << std::endl;
//                    if ( arms.top().id == 7374){
//                        break;
//                    }
//                    arms.pop();
//                }
//
//                std::cout << std::endl;
//                std::cout << "7374 = " << armsContainer[7374].lowerConfidenceBound << std::endl ;
//                break;
//            }

            bool bestArmFound;
            bestArmFound = iterationOfUCB();
            if (bestArmFound){
                topKArms.push_back(arms.top());
                arms.pop();
                bestArmCount++;
                std::cout << " Best arm number " << bestArmCount << " Position " << i <<std::endl;

                if (bestArmCount==numberOfBestArms)
                    break;
            }
        }
        if (bestArmCount!=numberOfBestArms){
            std::cout<< "UCB Stopped before reaching optimal" << std::endl;
        }
        std::cout << " Best arm number " << bestArmCount << " Position" << i<< " Max iter" << maxIterations << std::endl;
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
            //Checking if UCB should stop
            arms.push(bestArm);
//#ifdef DEBUG
            std::cout << "stopping UCB "<< std::setprecision (15)<< UCBofBestArm << "id " << bestArm.id <<  std::endl;
            std::cout << "stopping LCB "<< std::setprecision (15) <<  LCBofSecondBestArm << "id " << secondBestArm.id << std::endl;
            std::cout << "best "<< std::setprecision (15) <<  armsContainer[7374].lowerConfidenceBound
                    << " " << armsContainer[7374].estimateOfMean
                    << " " << armsContainer[7374].upperConfidenceBound
                    << " " << armsContainer[7374].numberOfPulls
                    << " p size" << arms.size()
                    << std::endl;

//#endif
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

    templateArm bestArm(){
        unsigned long bIndex = 0;
        float minTrueMean = armsContainer[0].trueMean();
        for (unsigned long i(1); i< armsContainer.size(); i++){
            float tmpTrueMean = armsContainer[i].trueMean();
            if ( tmpTrueMean < minTrueMean){
                minTrueMean = tmpTrueMean;
                bIndex = i;
            }
        }
        return armsContainer[bIndex];
    }
};


#endif //COMBINATORIAL_MAB_UCB_H
