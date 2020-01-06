#include "mchost.h"
#include "debug.h"


MCHost::MCHost(std::shared_ptr<ParameterStorage> parameterStorage) : parameterStorage(parameterStorage){
    DEBUG_FUNC_START

    hoppingSiteNumber=parameterStorage->parameters.at("hoppingSiteNumber");
    locLenA=parameterStorage->parameters.at("a");
    electrodeNumber=int(parameterStorage->electrodes.size());

    rates = new double*[hoppingSiteNumber];
    for(int i=0; i<hoppingSiteNumber;i++){
        rates[i]= new double[hoppingSiteNumber];
    }




    std::string dataFileName = parameterStorage->workingDirecotry+ "data.hdf5";
    dataFile = std::shared_ptr<DataFile>(new DataFile(dataFileName));

    voltageScanPointsNumber = (parameterStorage->parameters.at("voltageScanMax")-parameterStorage->parameters.at("voltageScanMin"))/parameterStorage->parameters.at("voltageScanResolution")+1;

    outputCurrentBuffer       = new double[voltageScanPointsNumber*voltageScanPointsNumber];
    outputCurrentUncertBuffer = new double[voltageScanPointsNumber*voltageScanPointsNumber];

    dataFile->createDataset("outputCurrent",       {voltageScanPointsNumber*voltageScanPointsNumber} );
    dataFile->createDataset("outputCurrentUncert", {voltageScanPointsNumber*voltageScanPointsNumber} );
    dataFile->createDataset("fitness",             {1});
    dataFile->createDataset("fitnessUncert",       {1});
    dataFile->createDataset("voltages",            {electrodeNumber} );

    
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


            // std::cout<<system->deltaEnergies[i][j]<<" ";
        }
        // std::cout<<std::endl;
    }
    // std::cout<<std::endl;

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


void MCHost::singleRun(){
    //run system until currents are in equilibrium
    DEBUG_FUNC_START

    // reset currents
    for (int i = 0; i < hoppingSiteNumber; i++){
        system->hoppingSites[i]->currentCounter    = 0;
        system->hoppingSites[i]->absCurrentCounter = 0;
    }
    // run equil steps
    int N=parameterStorage->parameters.at("equilSteps");
    for(int i=0; i<N;i++){
        system->calcEnergies();
        calcRates();
        makeSwap();
    }


    //split up in multiple runs to calc uncertainty of currents
    int uncertaintySplits=10;
    N=int(parameterStorage->parameters.at("calcCurrentSteps")/uncertaintySplits);

    outputCurrent     = 0;
    outputCurrentSqrt = 0;

    for(int j=0; j < uncertaintySplits; j++){
        // reset currents
        for (int i = 0; i < hoppingSiteNumber; i++){
            system->hoppingSites[i]->currentCounter    = 0;
            system->hoppingSites[i]->absCurrentCounter = 0;
        }
        // run pruductions steps
        for(int i=0; i<N;i++){
            system->calcEnergies();
            calcRates();
            makeSwap();
        }
        outputCurrent     +=           system->hoppingSites[parameterStorage->parameters.at("outputElectrode")+parameterStorage->parameters["acceptorNumber"]]->currentCounter;
        outputCurrentSqrt += std::pow(system->hoppingSites[parameterStorage->parameters.at("outputElectrode")+parameterStorage->parameters["acceptorNumber"]]->currentCounter,2);
    }

    outputCurrentStd=std::sqrt((outputCurrentSqrt-outputCurrent*outputCurrent/uncertaintySplits)/uncertaintySplits);
    outputCurrent/=uncertaintySplits;

    // //print currents
    // for (int i = 0; i < hoppingSiteNumber; i++){
    //     std::cout<<i<<" "<<system->hoppingSites[i]->absCurrentCounter<<std::endl;
    // }
   
    DEBUG_FUNC_END
}


void MCHost::runVoltageSetup(){
    DEBUG_FUNC_START

    for(int i=0; i < voltageScanPointsNumber; i++){            
        for(int j=0; j < voltageScanPointsNumber; j++){
            std::cout<<"scanning... "<<parameterStorage->parameters.at("voltageScanMin")+parameterStorage->parameters.at("voltageScanResolution")*i<<" "
                                     <<parameterStorage->parameters.at("voltageScanMin")+parameterStorage->parameters.at("voltageScanResolution")*j<<" ";
            system->setElectrodeVoltage(parameterStorage->parameters.at("inputElectrode1"),parameterStorage->parameters.at("voltageScanMin")+parameterStorage->parameters.at("voltageScanResolution")*i);
            system->setElectrodeVoltage(parameterStorage->parameters.at("inputElectrode2"),parameterStorage->parameters.at("voltageScanMin")+parameterStorage->parameters.at("voltageScanResolution")*j);

            system->updatePotential();

            singleRun();

            outputCurrentBuffer      [i*voltageScanPointsNumber+j]=outputCurrent;
            outputCurrentUncertBuffer[i*voltageScanPointsNumber+j]=outputCurrentStd;
            std:cout<<"current: "<<outputCurrentBuffer[i*voltageScanPointsNumber+j]<<std::endl;
        }
    }


    DEBUG_FUNC_END
}

void MCHost::saveResults(){
    DEBUG_FUNC_START

    double voltageBuffer[electrodeNumber]; 
    for(int i=0;i<electrodeNumber;i++){
        voltageBuffer[i]=parameterStorage->electrodes[i].voltage;
    }
    
    dataFile->addData("outputCurrent"      ,outputCurrentBuffer);
    dataFile->addData("outputCurrentUncert",outputCurrentUncertBuffer);
    dataFile->addData("fitness"            ,& fitness);
    dataFile->addData("fitnessUncert"      ,& fitnessUncert);
    dataFile->addData("voltages"           ,voltageBuffer);

    DEBUG_FUNC_END
}


void MCHost::calcFitness(){
    // el of [0,1]
    DEBUG_FUNC_START

    int maxIndex=0,minIndex=0;

    for(int i=0; i < voltageScanPointsNumber; i++){            
        for(int j=0; j < voltageScanPointsNumber; j++){
            if (outputCurrentBuffer[i*voltageScanPointsNumber+j]<outputCurrentBuffer[minIndex]){minIndex=i*voltageScanPointsNumber+j;}
            if (outputCurrentBuffer[i*voltageScanPointsNumber+j]>outputCurrentBuffer[maxIndex]){maxIndex=i*voltageScanPointsNumber+j;}
        }
    }

    fitness=0;
    fitnessUncert=0;
    double normed,desiredVal,normedUncert;
    for(int i=0; i < voltageScanPointsNumber; i++){            
        for(int j=0; j < voltageScanPointsNumber; j++){
            normed=(outputCurrentBuffer[i*voltageScanPointsNumber+j]-outputCurrentBuffer[minIndex])/(outputCurrentBuffer[maxIndex]-outputCurrentBuffer[minIndex]);
            normedUncert=std::sqrt(std::pow(outputCurrentUncertBuffer[i*voltageScanPointsNumber+j]/(outputCurrentBuffer[maxIndex]-outputCurrentBuffer[minIndex]),2)
                                  +std::pow((outputCurrentBuffer[i*voltageScanPointsNumber+j]-outputCurrentBuffer[minIndex])/std::pow(outputCurrentBuffer[maxIndex]-outputCurrentBuffer[minIndex],2)*outputCurrentUncertBuffer[maxIndex],2)
                                  +std::pow(((outputCurrentBuffer[i*voltageScanPointsNumber+j]-outputCurrentBuffer[minIndex])/std::pow(outputCurrentBuffer[maxIndex]-outputCurrentBuffer[minIndex],2)-1/(outputCurrentBuffer[maxIndex]-outputCurrentBuffer[minIndex]))*outputCurrentUncertBuffer[minIndex],2));
            desiredVal = desiredLogicFunction(parameterStorage->parameters.at("voltageScanMin")+parameterStorage->parameters.at("voltageScanResolution")*i,parameterStorage->parameters.at("voltageScanMin")+parameterStorage->parameters.at("voltageScanResolution")*j,parameterStorage->gate);
            fitness+=std::abs(normed-desiredVal);
            fitnessUncert+=normedUncert*normedUncert;
        }
    }
    fitness/=voltageScanPointsNumber*voltageScanPointsNumber;
    fitness=1-fitness;
    fitnessUncert=std::sqrt(fitnessUncert);
    fitnessUncert/=voltageScanPointsNumber*voltageScanPointsNumber;


    DEBUG_FUNC_END
}

void MCHost::optimizeMC(){
    DEBUG_FUNC_START
    std::cout<<"running optimization"<<std::endl;

    system->setElectrodeVoltage(parameterStorage->parameters.at("outputElectrode"),0);
    
    // init random voltages
    // for(int i=0;i<electrodeNumber;i++){
    //     if((i !=parameterStorage->parameters.at("outputElectrode")) & (i !=parameterStorage->parameters.at("inputElectrode1")) &(i !=parameterStorage->parameters.at("inputElectrode2"))){
    //         system->setElectrodeVoltage(i,enhance::random_double(parameterStorage->parameters.at("controlVoltageMin"),parameterStorage->parameters.at("controlVoltageMax")));
    //     }
    // }


    runVoltageSetup();
    calcFitness();
    saveResults();

    std::cout<<"fitness "<<fitness<<" +- "<<fitnessUncert<<std::endl;;

    double lastFitness=fitness;
    double lastVoltages[electrodeNumber];
    for(int i=0;i<electrodeNumber;i++){
        lastVoltages[i]=parameterStorage->electrodes[i].voltage;
    }


    for (size_t i = 0; i < 1000; i++){
        //get new random voltages
        std::cout<<"new random voltages: "<<std::endl;
        for(int i=0;i<electrodeNumber;i++){
            if((i !=parameterStorage->parameters.at("outputElectrode")) & (i !=parameterStorage->parameters.at("inputElectrode1")) &(i !=parameterStorage->parameters.at("inputElectrode2"))){
                system->setElectrodeVoltage(i,enhance::random_double(std::max(parameterStorage->parameters.at("controlVoltageMin"),parameterStorage->electrodes[i].voltage-parameterStorage->parameters.at("maxDeltaV")),std::min(parameterStorage->parameters.at("controlVoltageMax"),parameterStorage->electrodes[i].voltage+parameterStorage->parameters.at("maxDeltaV"))));
                std::cout<<i<<" "<<parameterStorage->electrodes[i].voltage<<std::endl;
            }
        }

        runVoltageSetup();
        calcFitness();
        saveResults();

        std::cout<<"fitness ("<<fitness<<" +- "<<fitnessUncert<<") lastFitness "<<lastFitness;
        if((fitness < lastFitness) & (enhance::fastExp((fitness-lastFitness)/parameterStorage->parameters.at("acceptanceFactor"))<enhance::random_double(0,1))){
            std::cout<<" not accepted "<<std::endl;
            //swap back
            for(int i=0;i<electrodeNumber;i++){
                system->setElectrodeVoltage(i,lastVoltages[i]);
            }
        }
        else{
            std::cout<<" accepted"<<std::endl;
            //setup for next iteration
            for(int i=0;i<electrodeNumber;i++){
                lastVoltages[i]=parameterStorage->electrodes[i].voltage;
            }
            lastFitness=fitness;
        }
    }
    

    DEBUG_FUNC_END
}


void MCHost::run(){
    DEBUG_FUNC_START
    std::cout<<"running fixed setup"<<std::endl;

    system->setElectrodeVoltage(parameterStorage->parameters.at("outputElectrode"),0);
    
    runVoltageSetup();
    calcFitness();
    saveResults();

    std::cout<<"fitness "<<fitness<<" +- "<<fitnessUncert<<std::endl;;


    DEBUG_FUNC_END
}



bool MCHost::desiredLogicFunction(double val1, double val2, std::string gate){
    DEBUG_FUNC_START

    bool b1 = val1 > parameterStorage->parameters.at("seperationVoltage");
    bool b2 = val2 > parameterStorage->parameters.at("seperationVoltage");

    if     (gate == "AND" ){ return  (b1 & b2);}
    else if(gate == "NAND"){ return !(b1 & b2);}
    else if(gate == "OR"  ){ return  (b1 | b2);}
    else if(gate == "NOR" ){ return !(b1 | b2);}
    else if(gate == "XOR" ){ return  (b1 ^ b2);}
    else if(gate == "NXOR"){ return !(b1 ^ b2);}
    else{
        throw std::runtime_error("logic operation not found");
    }

    DEBUG_FUNC_END
}