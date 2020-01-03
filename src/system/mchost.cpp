#include "mchost.h"
#include "debug.h"


MCHost::MCHost(std::shared_ptr<ParameterStorage> parameterStorage) : parameterStorage(parameterStorage){
    DEBUG_FUNC_START

    hoppingSiteNumber=this->parameterStorage->parameters.at("hoppingSiteNumber");
    locLenA=this->parameterStorage->parameters.at("a");
    currentConverganceIntervalCheckSteps=this->parameterStorage->parameters.at("currentConverganceIntervalCheckSteps");

    rates = new double*[hoppingSiteNumber];
    for(int i=0; i<hoppingSiteNumber;i++){
        rates[i]= new double[hoppingSiteNumber];
    }
    lastCurrentCounter    = new double[hoppingSiteNumber];
    lastAbsCurrentCounter = new double[hoppingSiteNumber];

    std::string dataFileName = parameterStorage->workingDirecotry+ "data.hdf5";
    dataFile = std::shared_ptr<DataFile>(new DataFile(dataFileName));

    dataFile->createDataset("currents",{hoppingSiteNumber});
    dataFile->createDataset("absCurrents",{hoppingSiteNumber});

    
    DEBUG_FUNC_END
}

void MCHost::setup(bool makeNewDevice)
{
    DEBUG_FUNC_START

    system.reset(new System(parameterStorage));
    system->initilizeMatrices();

    if(!makeNewDevice){
        system->loadDevice();
    }
    else{
        system->createRandomNewDevice();
    }
    
    system->getReadyForRun();

    DEBUG_FUNC_END
}

void MCHost::calcRates(){
    DEBUG_FUNC_START
    ratesSum=0;
    for(int i=0; i<hoppingSiteNumber;i++){
        for(int j=0; j<hoppingSiteNumber;j++){
            if(system->deltaEnergies[i][j] < 0){ 
                rates[i][j]=enhance::mediumFastExp(-2*system->distances[i][j]/locLenA);
                // std::cout<<rates[i][j]<<" < 0 "<<i<<" "<<j<<" r/a "<<system->distances[i][j]/locLenA<<std::endl;
            }
            else if(system->deltaEnergies[i][j] >0){
                rates[i][j]=enhance::mediumFastExp(-2*system->distances[i][j]/locLenA-system->deltaEnergies[i][j]);
                // std::cout<<rates[i][j]<<" > 0 "<<i<<" "<<j<<" r/a "<<-2*system->distances[i][j]/locLenA<<" dE "<<system->deltaEnergies[i][j]<<std::endl;
            }
            else{ // !!!!!!!!!! BUG POTENTIAL if deltaE==0 by accident -> rate=0, but should be 1 !!!!!!!!!!
                rates[i][j]=0;
                // std::cout<<rates[i][j]<<" = 0 "<<i<<" "<<j<<std::endl;
            }
            if (std::isinf(rates[i][j])){
                // std::cout<<"rate inf "<<rates[i][j]<<" i "<<i<<" j "<<j<<std::endl;
            }
            ratesSum+=rates[i][j];
        }
    }    
    DEBUG_FUNC_END
}


void MCHost::makeSwap(){
    DEBUG_FUNC_START
    double rndNumber=enhance::random_double(0,ratesSum);
    double partRatesSum=0;
    for(int i=0; i<hoppingSiteNumber;i++){
        for(int j=0; j<hoppingSiteNumber;j++){
            partRatesSum+=rates[i][j];
            // std::cout<<partRatesSum<<" "<<rates[i][j]<<" "<<rndNumber<<" "<<ratesSum<<std::endl;
            if(partRatesSum > rndNumber){
                system->hoppingSites[i]->setOccupation(false);
                system->hoppingSites[i]->currentCounter++;
                system->hoppingSites[i]->absCurrentCounter++;
                system->hoppingSites[j]->setOccupation(true);
                system->hoppingSites[j]->currentCounter--;
                system->hoppingSites[j]->absCurrentCounter++;

                // std::cout<<"swapped "<<i<<" "<<j<<" "<<setw(9);
                goto endDoubleLoop;
            }
        }
        if(i== hoppingSiteNumber-1){
            // std::cout<<"no swapp found!"<<partRatesSum<<" "<<rndNumber<<" "<<ratesSum<<" ";

        }
    }    
    endDoubleLoop:;
    DEBUG_FUNC_END
}


void MCHost::singleRun(int N){
    //run system until currents are in equilibrium
    DEBUG_FUNC_START

    for(int i=0; i<N;i++){
        system->calcEnergies();
        calcRates();
        makeSwap();

        // // print occupations
        // std::cout<<i<<" ";
        // for (int i = 0; i < hoppingSiteNumber; i++){
        //     std::cout<<setw(2)<<" "<<system->hoppingSites[i]->getOccupation();
        // }
        // std::cout<<std::endl;

    
        if((i-1)%currentConverganceIntervalCheckSteps==0){

            // print currents
            std::cout<<"currents "<<i<<std::endl;
            for (int i = 0; i < hoppingSiteNumber; i++){
                std::cout<<i<<" "<<system->hoppingSites[i]->currentCounter<<" "<<system->hoppingSites[i]->absCurrentCounter<<std::endl;
 
            }
            std::cout<<std::endl;



            for (int i = 0; i < hoppingSiteNumber; i++){
                lastCurrentCounter[i]   =system->hoppingSites[i]->currentCounter;
                lastAbsCurrentCounter[i]=system->hoppingSites[i]->absCurrentCounter;

                system->hoppingSites[i]->currentCounter    = 0;
                system->hoppingSites[i]->absCurrentCounter = 0;
            }

            dataFile->addData("currents"   ,lastCurrentCounter);
            dataFile->addData("absCurrents",lastAbsCurrentCounter);
        }
    }

    DEBUG_FUNC_END
}


void MCHost::run(){
    DEBUG_FUNC_START

    singleRun(parameterStorage->parameters.at("steps"));

    // system->setElectrodeVoltage(0,0);
    // system->setElectrodeVoltage(1,50);
    system->updatePotential();

    // for (int i = 0; i < hoppingSiteNumber; i++){
    //     system->hoppingSites[i]->currentCounter=0;
    //     system->hoppingSites[i]->absCurrentCounter=0;
    // }

    // singleRun(parameterStorage->parameters.at("steps"));


    DEBUG_FUNC_END
}

// void MCHost::printCurrents(int step){


// }
