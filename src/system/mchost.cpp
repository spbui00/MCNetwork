#include "mchost.h"
#include "debug.h"


MCHost::MCHost(std::shared_ptr<ParameterStorage> parameterStorage) : parameterStorage(parameterStorage){
    DEBUG_FUNC_START

    hoppingSiteNumber=parameterStorage->parameters.at("hoppingSiteNumber");
    electrodeNumber=int(parameterStorage->electrodes.size());


    std::string dataFileName = parameterStorage->workingDirecotry+ "data.hdf5";
    dataFile = std::shared_ptr<DataFile>(new DataFile(dataFileName));

    voltageScanPoints = parameterStorage->parameters.at("voltageScanPoints");
    voltageScanResolution=(parameterStorage->parameters.at("voltageScanMax")-parameterStorage->parameters.at("voltageScanMin"))/double(voltageScanPoints-1);

    outputCurrentBuffer       = new double[voltageScanPoints*voltageScanPoints];
    outputCurrentUncertBuffer = new double[voltageScanPoints*voltageScanPoints];

    dataFile->createDataset("outputCurrent",       {voltageScanPoints,voltageScanPoints} );
    dataFile->createDataset("outputCurrentUncert", {voltageScanPoints,voltageScanPoints} );
    dataFile->createDataset("fitness",             {1});
    dataFile->createDataset("fitnessUncert",       {1});
    dataFile->createDataset("optEnergy",           {1});
    dataFile->createDataset("voltages",            {electrodeNumber} );

    
    DEBUG_FUNC_END
}

void MCHost::setup(bool makeNewDevice)
{
    DEBUG_FUNC_START

    // systems.push_back(std::unique_ptr<System>(new System(parameterStorage)));
    systems.push_back( new System(parameterStorage));
    systems[0]->initilizeMatrices();

    if(!makeNewDevice){
        systems[0]->loadDevice();
    }
    else{
        systems[0]->createRandomNewDevice();
    }
    
    systems[0]->getReadyForRun();

    for(int i=1; i < parameterStorage->parameters.at("threads"); i++){
        systems.push_back(new System(*systems[0],bool(parameterStorage->parameters.at("shareMemory"))));
    }

    DEBUG_FUNC_END
}



void MCHost::singleRun(){
    //run system until currents are in equilibrium
    DEBUG_FUNC_START

    // *(systems[0]->storeKnownStates)=true;
    std::vector<std::thread> threads;

    bool multiProcessingMode = parameterStorage->parameters.at("threads") > 1 ;

    // run equil steps
    int N=parameterStorage->parameters.at("equilSteps");
    if (multiProcessingMode){
        for(int k=0; k < parameterStorage->parameters.at("threads"); k++){
            threads.push_back(std::thread(&System::run, std::ref(*systems[k]),N));
        }
        for(int k=0; k < parameterStorage->parameters.at("threads"); k++){
            threads[k].join();
        }
        threads.clear();
    }
    else{
        systems[0]->run(N);
    }



    //split up in multiple runs to calc uncertainty of currents
    int runsPerThread=std::ceil(20.0/parameterStorage->parameters.at("threads"));

    N=int(parameterStorage->parameters.at("calcCurrentSteps")/(runsPerThread*parameterStorage->parameters.at("threads")));

    outputCurrent     = 0;
    outputCurrentSqrt = 0;



    for(int j=0; j < runsPerThread; j++){
        for(int k=0; k < parameterStorage->parameters.at("threads"); k++){
            // reset currents
            systems[k]->time=0;
            for (int i = 0; i < hoppingSiteNumber; i++){
                // std::cout<<i<<" "<<systems[k]->currentCounter[k]<<" "<<systems[k]->time<<std::endl;
                systems[k]->currentCounter[i] = 0;
            }
            
            // run pruductions steps
            if (multiProcessingMode){
                threads.push_back(std::thread(&System::run, systems[k],N));
            }
            else{
                systems[k]->run(N); // anyway k == 0
            }


        }

        for(int k=0; k < parameterStorage->parameters.at("threads"); k++){
            if (multiProcessingMode){
                threads[k].join();
            }

            outputCurrent     +=          systems[k]->currentCounter[int(parameterStorage->parameters.at("outputElectrode")+parameterStorage->parameters["acceptorNumber"])]/systems[k]->time;
            outputCurrentSqrt += std::pow(systems[k]->currentCounter[int(parameterStorage->parameters.at("outputElectrode")+parameterStorage->parameters["acceptorNumber"])]/systems[k]->time,2);
            // std::cout<<"curr: "<<" "<<outputCurrent/(j+1) <<" +- "<<std::sqrt((outputCurrentSqrt-outputCurrent*outputCurrent/(j+1)))/(j+1) <<std::endl;
            // std::cout<<"thread: "<<k<<" run: "<<j<<" curr: "<<" "<<systems[k]->currentCounter[int(parameterStorage->parameters.at("outputElectrode")+parameterStorage->parameters["acceptorNumber"])]/systems[k]->time <<std::endl;
        }
        if (multiProcessingMode){
            threads.clear();
        }
    }


    outputCurrentStd=std::sqrt((outputCurrentSqrt-outputCurrent*outputCurrent/(runsPerThread*parameterStorage->parameters.at("threads"))))/(runsPerThread*parameterStorage->parameters.at("threads"));
    outputCurrent/=(runsPerThread*parameterStorage->parameters.at("threads"));

    // //print currents
    // for (int i = 0; i < hoppingSiteNumber; i++){
    //     std::cout<<i<<" "<<system->hoppingSites[i]->absCurrentCounter<<std::endl;
    // }
   
    DEBUG_FUNC_END
}



void MCHost::runVoltageSetup(){
    DEBUG_FUNC_START

    for(int i=0; i < voltageScanPoints; i++){            
        for(int j=0; j < voltageScanPoints; j++){
            std::cout<<"scanning... "<<parameterStorage->parameters.at("voltageScanMin")+voltageScanResolution*i<<" "
                                     <<parameterStorage->parameters.at("voltageScanMin")+voltageScanResolution*j<<" ";
            parameterStorage->electrodes[parameterStorage->parameters.at("inputElectrode1")].voltage=parameterStorage->parameters.at("voltageScanMin")+voltageScanResolution*i;
            parameterStorage->electrodes[parameterStorage->parameters.at("inputElectrode2")].voltage=parameterStorage->parameters.at("voltageScanMin")+voltageScanResolution*j;


            //recalc potential
            for(int k=0; k < parameterStorage->parameters.at("threads"); k++){
                systems[k]->resetPotential();
            }
            systems[0]->recalcPotential();
            for(int k=0; k < parameterStorage->parameters.at("threads"); k++){
                systems[k]->setNewPotential();
            }

            //reset stored states
            if (parameterStorage->parameters.at("shareMemory")){
                std::cout<<"maximal size of stored states: "<<(hoppingSiteNumber*hoppingSiteNumber+1)*8/1e6*systems[0]->knownRatesSum->size()<<" mb; "<<systems[0]->knownRatesSum->size()<<" states"<<std::endl;
                systems[0]->konwnPartRatesSumList->clear();
                systems[0]->knownRatesSum        ->clear();
            }
            else{
                for(int k=0; k < parameterStorage->parameters.at("threads"); k++){
                    std::cout<<"maximal size of stored states: "<<(hoppingSiteNumber*hoppingSiteNumber+1)*8/1e6*systems[k]->knownRatesSum->size()<<" mb; "<<systems[k]->knownRatesSum->size()<<" states"<<std::endl;
                    systems[k]->konwnPartRatesSumList->clear();
                    systems[k]->knownRatesSum        ->clear();
                }
            }

            singleRun();

            outputCurrentBuffer      [i*voltageScanPoints+j]=outputCurrent;
            outputCurrentUncertBuffer[i*voltageScanPoints+j]=outputCurrentStd;
            std:cout<<"current: "<<outputCurrentBuffer[i*voltageScanPoints+j]<<" +- "<<outputCurrentUncertBuffer[i*voltageScanPoints+j]<<std::endl;
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
    dataFile->addData("optEnergy"          ,& optEnergy);
    dataFile->addData("voltages"           ,voltageBuffer);

    DEBUG_FUNC_END
}


void MCHost::calcOptimizationEnergy(){
    // el of [0,1]
    DEBUG_FUNC_START

    int maxIndex=0,minIndex=0;

    for(int i=0; i < voltageScanPoints; i++){            
        for(int j=0; j < voltageScanPoints; j++){
            if (outputCurrentBuffer[i*voltageScanPoints+j]<outputCurrentBuffer[minIndex]){minIndex=i*voltageScanPoints+j;}
            if (outputCurrentBuffer[i*voltageScanPoints+j]>outputCurrentBuffer[maxIndex]){maxIndex=i*voltageScanPoints+j;}
        }
    }

    fitness=0;
    fitnessUncert=0;
    double normed,desiredVal,normedUncert;
    for(int i=0; i < voltageScanPoints; i++){            
        for(int j=0; j < voltageScanPoints; j++){
            normed=(outputCurrentBuffer[i*voltageScanPoints+j]-outputCurrentBuffer[minIndex])/(outputCurrentBuffer[maxIndex]-outputCurrentBuffer[minIndex]);
            normedUncert=std::sqrt(std::pow(outputCurrentUncertBuffer[i*voltageScanPoints+j]/(outputCurrentBuffer[maxIndex]-outputCurrentBuffer[minIndex]),2)
                                  +std::pow((outputCurrentBuffer[i*voltageScanPoints+j]-outputCurrentBuffer[minIndex])/std::pow(outputCurrentBuffer[maxIndex]-outputCurrentBuffer[minIndex],2)*outputCurrentUncertBuffer[maxIndex],2)
                                  +std::pow(((outputCurrentBuffer[i*voltageScanPoints+j]-outputCurrentBuffer[minIndex])/std::pow(outputCurrentBuffer[maxIndex]-outputCurrentBuffer[minIndex],2)-1/(outputCurrentBuffer[maxIndex]-outputCurrentBuffer[minIndex]))*outputCurrentUncertBuffer[minIndex],2));
            desiredVal = desiredLogicFunction(parameterStorage->parameters.at("voltageScanMin")+voltageScanResolution*i,parameterStorage->parameters.at("voltageScanMin")+voltageScanResolution*j,parameterStorage->gate);
            fitness+=std::abs(normed-desiredVal);
            fitnessUncert+=normedUncert*normedUncert;
        }
    }
    fitness/=voltageScanPoints*voltageScanPoints;
    fitness=1-fitness;
    fitnessUncert=std::sqrt(fitnessUncert);
    fitnessUncert/=voltageScanPoints*voltageScanPoints;
    normedDiff= (outputCurrentBuffer[maxIndex]-outputCurrentBuffer[minIndex])/(2*std::max(std::abs(outputCurrentBuffer[maxIndex]),std::abs(outputCurrentBuffer[minIndex])));

    optEnergy = fitness - fitnessUncert*parameterStorage->parameters.at("fitnessUncertWeight") + normedDiff*parameterStorage->parameters.at("diffWeight");


    DEBUG_FUNC_END
}

void MCHost::optimizeMC(bool rndStart /*= false*/){
    DEBUG_FUNC_START
    std::cout<<"running optimization"<<std::endl;

    parameterStorage->electrodes[parameterStorage->parameters.at("outputElectrode")].voltage=0;
    
    // init random voltages
    if (rndStart){
        std::cout<<"new random voltages: "<<std::endl;
        for(int i=0;i<electrodeNumber;i++){
            if((i !=parameterStorage->parameters.at("outputElectrode")) & (i !=parameterStorage->parameters.at("inputElectrode1")) &(i !=parameterStorage->parameters.at("inputElectrode2"))){
                parameterStorage->electrodes[i].voltage=enhance::random_double(parameterStorage->parameters.at("controlVoltageMin"),parameterStorage->parameters.at("controlVoltageMax"));
                std::cout<<i<<" "<<parameterStorage->electrodes[i].voltage<<std::endl;
            }
        }
    }

    auto startTime = chrono::steady_clock::now();

    runVoltageSetup();
    calcOptimizationEnergy();
    saveResults();

    std::cout<<"optEnergy: "<<optEnergy    <<" fitness: ("<<fitness    <<" +- "<<fitnessUncert    <<") normedDiff: "<<normedDiff<<std::endl;;
    auto endTime = chrono::steady_clock::now();
    std::cout << "time per VoltageSetup = " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count()/1000.0 << " s" << std::endl;
    



    double lastFitness       = fitness;
    double lastFitnessUncert = fitnessUncert;
    double lastOptEnergy     = optEnergy;
    double lastNormedDiff    = normedDiff;
    double lastVoltages[electrodeNumber];
    for(int i=0;i<electrodeNumber;i++){
        lastVoltages[i]=parameterStorage->electrodes[i].voltage;
    }

    int maxIncreases=4;
    int increaseNumber=0;
    for (size_t i = 0; i < 1000000; i++){
        auto startTime = chrono::steady_clock::now();

        //get new random voltages
        std::cout<<"new random voltages: "<<std::endl;
        for(int i=0;i<electrodeNumber;i++){
            if((i !=parameterStorage->parameters.at("outputElectrode")) & (i !=parameterStorage->parameters.at("inputElectrode1")) &(i !=parameterStorage->parameters.at("inputElectrode2"))){
                parameterStorage->electrodes[i].voltage=enhance::random_double(std::max(parameterStorage->parameters.at("controlVoltageMin"),parameterStorage->electrodes[i].voltage-parameterStorage->parameters.at("maxDeltaV")),std::min(parameterStorage->parameters.at("controlVoltageMax"),parameterStorage->electrodes[i].voltage+parameterStorage->parameters.at("maxDeltaV")));
                std::cout<<i<<" "<<parameterStorage->electrodes[i].voltage<<std::endl;
            }
        }

        runVoltageSetup();
        calcOptimizationEnergy();
        saveResults();

        std::cout<<"now:  optEnergy: "<<optEnergy    <<" fitness: ("<<fitness    <<" +- "<<fitnessUncert    <<") normedDiff: "<<normedDiff<<"\nlast: optEnergy: "<<lastOptEnergy<<" fitness: ("<<lastFitness<<" +- "<<lastFitnessUncert<<") normedDiff: "<<lastNormedDiff<<std::endl;
        if((optEnergy < lastOptEnergy) & (enhance::fastExp((optEnergy-lastOptEnergy)/parameterStorage->parameters.at("acceptanceFactor"))<enhance::random_double(0,1))){
            std::cout<<"-- not accepted --"<<std::endl;
            //swap back
            for(int i=0;i<electrodeNumber;i++){
                parameterStorage->electrodes[i].voltage=lastVoltages[i];
            }
        }
        else{
            std::cout<<"-- accepted --"<<std::endl;
            //setup for next iteration
            for(int i=0;i<electrodeNumber;i++){
                lastVoltages[i]=parameterStorage->electrodes[i].voltage;
            }
            lastFitness       = fitness;
            lastFitnessUncert = fitnessUncert;
            lastOptEnergy     = optEnergy;
            lastNormedDiff    = normedDiff;

        }
        auto endTime = chrono::steady_clock::now();
        std::cout << "time per VoltageSetup = " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count()/1000.0 << " s" << std::endl;
    

        if((increaseNumber < maxIncreases) and (fitness+fitnessUncert*2)>1){
            increaseNumber++;
            parameterStorage->parameters["calcCurrentSteps"]*=2;
            std::cout<<"############ steps increased!! now: "<<parameterStorage->parameters["calcCurrentSteps"]<<" #############"<<std::endl;
        }
    }
    


    DEBUG_FUNC_END
}


void MCHost::run(){
    DEBUG_FUNC_START
    std::cout<<"running fixed setup"<<std::endl;

    auto startTime = chrono::steady_clock::now();

    parameterStorage->electrodes[parameterStorage->parameters.at("outputElectrode")].voltage=0;

    runVoltageSetup();
    calcOptimizationEnergy();
    saveResults();

    std::cout<<"optEnergy: "<<optEnergy    <<" fitness: ("<<fitness    <<" +- "<<fitnessUncert    <<") normedDiff: "<<normedDiff<<std::endl;;

    auto endTime = chrono::steady_clock::now();
    std::cout << "time elapsed = " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count()/1000.0 << " s" << std::endl;
    
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