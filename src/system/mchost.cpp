#include "mchost.h"
#include "debug.h"


MCHost::MCHost(std::shared_ptr<ParameterStorage> parameterStorage) : parameterStorage(parameterStorage){
    DEBUG_FUNC_START

    hoppingSiteNumber=this->parameterStorage->parameters.at("hoppingSiteNumber");
    locLenA=this->parameterStorage->parameters.at("a");

    rates = new double*[hoppingSiteNumber];
    for(int i; i<hoppingSiteNumber;i++){
        rates[i]= new double[hoppingSiteNumber];
    }
    
    DEBUG_FUNC_END
}

void MCHost::setup(std::string deviceFileName /*=""*/)
{
    DEBUG_FUNC_START

    system.reset(new System(parameterStorage));
    system->initilizeMatrices();

    if(deviceFileName!=""){
        system->loadDevice(deviceFileName);
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
                rates[i][j]=enhance::fastExp(-2*system->distances[i][j]/locLenA);
                // std::cout<<rates[i][j]<<" < 0 "<<i<<" "<<j<<" r/a "<<system->distances[i][j]/locLenA<<std::endl;
            }
            else if(system->deltaEnergies[i][j] >0){
                rates[i][j]=enhance::fastExp(-2*system->distances[i][j]/locLenA-system->deltaEnergies[i][j]);
                // std::cout<<rates[i][j]<<" > 0 "<<i<<" "<<j<<" r/a "<<-2*system->distances[i][j]/locLenA<<" dE "<<system->deltaEnergies[i][j]<<std::endl;
            }
            else{ // !!!!!!!!!! BUG POTENTIAL if deltaE==0 by accident -> rate=0, but should be 1 !!!!!!!!!!
                rates[i][j]=0;
                // std::cout<<rates[i][j]<<" = 0 "<<i<<" "<<j<<std::endl;
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
                system->hoppingSites[j]->setOccupation(true);
                std::cout<<"swapped "<<i<<" "<<j<<" "<<setw(9);
                goto endDoubleLoop;
            }
        }
    }    
    endDoubleLoop:;
    DEBUG_FUNC_END
}


void MCHost::singleRun(int N){
    DEBUG_FUNC_START

    for(int i=0; i<N;i++){
        system->calcEnergies();
        calcRates();
        makeSwap();

        // std::cout<<i<<" ";
        // for (int i = 0; i < hoppingSiteNumber; i++){
        //     std::cout<<setw(2)<<" "<<system->hoppingSites[i]->getOccupation();
        // }
        // std::cout<<std::endl;
    }
    DEBUG_FUNC_END
}


void MCHost::run(){
    DEBUG_FUNC_START

    singleRun(parameterStorage->parameters.at("steps"));

    system->setElectrodeVoltage(0,0);
    system->setElectrodeVoltage(1,0);
    system->updatePotential();

    singleRun(parameterStorage->parameters.at("steps"));


    DEBUG_FUNC_END
}
