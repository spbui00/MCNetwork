#include "mchost.h"
#include "debug.h"


MCHost::MCHost(std::shared_ptr<ParameterStorage> _parameterStorage)
{
    DEBUG_FUNC_START

    parameterStorage=_parameterStorage;
    acceptorNumber=parameterStorage->parameters.at("acceptorNumber");
    locLenA=parameterStorage->parameters.at("a");

    rates = new double*[acceptorNumber];
    for(int i; i<acceptorNumber;i++){
        rates[i]= new double[acceptorNumber];
    }
    
    DEBUG_FUNC_END
}

void MCHost::setup()
{
    DEBUG_FUNC_START

    
    system.reset(new System(parameterStorage));
    system->setup();

    DEBUG_FUNC_END
}

void MCHost::calcRates(){
    DEBUG_FUNC_START
    ratesSum=0;
    for(int i=0; i<acceptorNumber;i++){
        for(int j=0; j<acceptorNumber;j++){
            if(system->deltaEnergies[i][j] < 0){ 
                rates[i][j]=enhance::fastExp(-2*system->distances[i][j]/locLenA);
            }
            else if(system->deltaEnergies[i][j] >0){
                rates[i][j]=enhance::fastExp(-2*system->distances[i][j]/locLenA-system->deltaEnergies[i][j]);
            }
            else{ // !!!!!!!!!! BUG POTENTIAL if deltaE==0 by accident -> rate=0, but should be 1 !!!!!!!!!!
                rates[i][j]=0;
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
    for(int i=0; i<acceptorNumber;i++){
        for(int j=0; j<acceptorNumber;j++){
            partRatesSum+=rates[i][j];
            if(partRatesSum > rndNumber){
                system->dopants[i].occupied=false;
                system->dopants[j].occupied=true;
                goto endDoubleLoop;
            }
        }
    }    
    endDoubleLoop:;
    DEBUG_FUNC_END
}


void MCHost::run(int N){
    DEBUG_FUNC_START
    for(int i=0; i<N;i++){
        system->calcEnergies();
        calcRates();
        makeSwap();
    }
    DEBUG_FUNC_END
}