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


    std::vector<templateArm> armsContainer;
    std::priority_queue<templateArm, std::vector<templateArm>, std::greater<templateArm> > arms;
    std::vector<unsigned long > topKArms;
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
//                std::this_thread::sleep_for(700);
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
//        std::cout << numberOfThreads << " number of threads used in each batch of initialization.\n";
        std::vector<std::thread> initThreads(numberOfThreads);

        unsigned long chunkSize = (numberOfArms/numberOfThreads);
        for(unsigned t = 0; t < numberOfThreads; t++){
            unsigned long armIndexStart = t*chunkSize;
            unsigned long armIndexEnd = (t+1)*chunkSize;

//            std::cout << "Starting thread for arm " << armIndexStart
//                      << " to " << armIndexEnd << std::endl;
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

    void storeTopKArms(){
        for (int i=0; i<20; i++) {
            topKArms.push_back(arms.top().id);
            arms.pop();
        }
    }
    void runUCB(unsigned long maxIterations){
        for (unsigned long i(0); i < maxIterations; i++){
            bool stop;
            stop = iterationOfUCB();
            if (stop){
                storeTopKArms();
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


#endif //COMBINATORIAL_MAB_UCB_H
