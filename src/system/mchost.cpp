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
                *(systems[0]->storeKnownStates)=true;
            }
            else{
                for(int k=0; k < parameterStorage->parameters.at("threads"); k++){
                    std::cout<<"maximal size of stored states: "<<(hoppingSiteNumber*hoppingSiteNumber+1)*8/1e6*systems[k]->knownRatesSum->size()<<" mb; "<<systems[k]->knownRatesSum->size()<<" states"<<std::endl;
                    systems[k]->konwnPartRatesSumList->clear();
                    systems[k]->knownRatesSum        ->clear();
                    *(systems[k]->storeKnownStates)=true;
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

    dataFile->createDataset("accepted",{1});
    double accepted=-1;

    std::cout<<"running MC optimization"<<std::endl;

    parameterStorage->electrodes[parameterStorage->parameters.at("outputElectrode")].voltage=0;
    
    // init random voltages
    if (rndStart){
        std::cout<<"------ searching for start point ------"<<std::endl;

        std::vector<int> controlElectrodeIndices;
        for(int i=0;i<electrodeNumber;i++){
            if((i !=parameterStorage->parameters.at("outputElectrode")) & (i !=parameterStorage->parameters.at("inputElectrode1")) &(i !=parameterStorage->parameters.at("inputElectrode2"))){
                controlElectrodeIndices.push_back(i);
            }
        }
        int controlElectrodes=controlElectrodeIndices.size();
        
        std::vector<std::pair<std::vector<double>,double>> voltagesSet; //pair of control voltages and optEnergy
        auto voltagesSetComparator = []( const std::pair<std::vector<double>,double>& l, const std::pair<std::vector<double>,double>& r) { return l.second > r.second; }; //needed to find best start point

        for(int k=0; k < parameterStorage->parameters.at("MCOptStartPoints"); k++){
            voltagesSet.push_back(std::pair<std::vector<double>,double>(std::vector<double>(controlElectrodes),0));
            for(int j=0; j < controlElectrodes; j++){
                voltagesSet[k].first[j]=enhance::random_double(parameterStorage->parameters.at("controlVoltageMin"),parameterStorage->parameters.at("controlVoltageMax"));
            }
        }

        //run rnd start candidates
        for(int k=0; k < parameterStorage->parameters.at("MCOptStartPoints"); k++){
            std::cout<<"rnd start point: "<<k<<" voltages:";
            for(int i=0; i < controlElectrodes; i++){
                parameterStorage->electrodes[controlElectrodeIndices[i]].voltage=voltagesSet[k].first[i];
                std::cout<<" "<<controlElectrodeIndices[i]<<": "<<parameterStorage->electrodes[controlElectrodeIndices[i]].voltage;
            }
            std::cout<<std::endl;

            auto startTime = chrono::steady_clock::now();

            runVoltageSetup();
            calcOptimizationEnergy();
            saveResults();
            dataFile->addData("accepted",& accepted);
            
            voltagesSet[k].second=optEnergy;


            std::cout<<"optEnergy: "<<optEnergy    <<" fitness: ("<<fitness    <<" +- "<<fitnessUncert    <<") normedDiff: "<<normedDiff<<std::endl;
            auto endTime = chrono::steady_clock::now();
            std::cout << "time per VoltageSetup = " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count()/1000.0 << " s" << std::endl;
        }

        std::cout<<"start search done! sorted results: "<<std::endl;
        std::sort(voltagesSet.begin(),voltagesSet.end(),voltagesSetComparator);
        for(int k=0; k < parameterStorage->parameters.at("MCOptStartPoints"); k++){
            std::cout<<"genome: "<<k+1<<" optEnergy: "<<voltagesSet[k].second<<" voltages: ";
            for(int i=0; i < controlElectrodes; i++){
                std::cout<<" "<<controlElectrodeIndices[i]<<": "<<voltagesSet[k].first[i];
            }
            std::cout<<std::endl;
        }

        //set best start
        for(int i=0; i < controlElectrodes; i++){
            parameterStorage->electrodes[controlElectrodeIndices[i]].voltage=voltagesSet[0].first[i];
        }
    }


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
            accepted=0;
            //swap back
            for(int i=0;i<electrodeNumber;i++){
                parameterStorage->electrodes[i].voltage=lastVoltages[i];
            }
        }
        else{
            std::cout<<"-- accepted --"<<std::endl;
            accepted=1;
            //setup for next iteration
            for(int i=0;i<electrodeNumber;i++){
                lastVoltages[i]=parameterStorage->electrodes[i].voltage;
            }
            lastFitness       = fitness;
            lastFitnessUncert = fitnessUncert;
            lastOptEnergy     = optEnergy;
            lastNormedDiff    = normedDiff;

        }
        dataFile->addData("accepted",& accepted);


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


void MCHost::optimizeGenetic(){
    DEBUG_FUNC_START

    dataFile->createDataset("generation",{1});

    // lambda function needed later to sort genome set
    auto genomeComparator = []( const std::pair<std::vector<double>,double>& l, const std::pair<std::vector<double>,double>& r) { return l.second > r.second; };


    std::cout<<"running optimization - genetic"<<std::endl;

    parameterStorage->electrodes[parameterStorage->parameters.at("outputElectrode")].voltage=0;

    std::vector<int> controlElectrodeIndices;
    for(int i=0;i<electrodeNumber;i++){
        if((i !=parameterStorage->parameters.at("outputElectrode")) & (i !=parameterStorage->parameters.at("inputElectrode1")) &(i !=parameterStorage->parameters.at("inputElectrode2"))){
            controlElectrodeIndices.push_back(i);
        }
    }
    int controlElectrodes=controlElectrodeIndices.size();

    std::vector<std::pair<std::vector<double>,double>> genomeSet; //pair of control voltages and optEnergy
    for(int i=0; i < 25; i++){
        genomeSet.push_back(std::pair<std::vector<double>,double>(std::vector<double>(controlElectrodes),0));
        for(int j=0; j < controlElectrodes; j++){
            genomeSet[i].first[j]=enhance::random_double(parameterStorage->parameters.at("controlVoltageMin"),parameterStorage->parameters.at("controlVoltageMax"));
        }

    }

    //run first generation
    double generation=1;
    std::cout<<"------------------------------ run geneartion "<<generation<< " ------------------------------"<<std::endl;
    for(int k=0; k < 25; k++){
        std::cout<<"genome: "<<k<<" voltages:";
        for(int i=0; i < controlElectrodes; i++){
            parameterStorage->electrodes[controlElectrodeIndices[i]].voltage=genomeSet[k].first[i];
            std::cout<<" "<<controlElectrodeIndices[i]<<": "<<parameterStorage->electrodes[controlElectrodeIndices[i]].voltage;
        }
        std::cout<<std::endl;

        auto startTime = chrono::steady_clock::now();

        runVoltageSetup();
        calcOptimizationEnergy();
        saveResults();
        dataFile->addData("generation",& generation);
        
        genomeSet[k].second=optEnergy;


        std::cout<<"optEnergy: "<<optEnergy    <<" fitness: ("<<fitness    <<" +- "<<fitnessUncert    <<") normedDiff: "<<normedDiff<<std::endl;
        auto endTime = chrono::steady_clock::now();
        std::cout << "time per VoltageSetup = " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count()/1000.0 << " s" << std::endl;
    }

    std::cout<<"generation "<<generation<< " done! sorted results: "<<std::endl;
    std::sort(genomeSet.begin(),genomeSet.end(),genomeComparator);
    for(int k=0; k < 25; k++){
        std::cout<<"genome: "<<k+1<<" optEnergy: "<<genomeSet[k].second<<" voltages: ";
        for(int i=0; i < controlElectrodes; i++){
            std::cout<<" "<<controlElectrodeIndices[i]<<": "<<genomeSet[k].first[i];
        }
        std::cout<<std::endl;
    }

    while(true){
        // ---------- setup next generation -------------
        generation++;

        //first 5 genomes dont need to be changed

        //genome 6-10, rnd Bias
        for(int k=5; k < 10; k++){
            for(int i=0; i < controlElectrodes; i++){
                genomeSet[k].first[i]=enhance::random_double(std::max(parameterStorage->parameters.at("controlVoltageMin"),genomeSet[k-5].first[i]-parameterStorage->parameters.at("maxDeltaV")),std::min(parameterStorage->parameters.at("controlVoltageMax"),genomeSet[k-5].first[i]+parameterStorage->parameters.at("maxDeltaV")));
            }
            genomeSet[k].second=0;
        }

        //genome 11-15, crossover
        for(int k=10; k < 15; k++){
            for(int i=0; i < controlElectrodes; i++){
                if (enhance::random_double(0,1)>0.5){
                    genomeSet[k].first[i]=genomeSet[k-10].first[i];
                }
                else{
                    genomeSet[k].first[i]=genomeSet[k-9].first[i];
                }
            }
            genomeSet[k].second=0;
        }

        //genome 16-20, rnd crossover
        for(int k=15; k < 20; k++){
            for(int i=0; i < controlElectrodes; i++){
                if (enhance::random_double(0,1)>0.5){
                    genomeSet[k].first[i]=genomeSet[k-15].first[i];
                }
                else{
                    genomeSet[k].first[i]=enhance::random_double(parameterStorage->parameters.at("controlVoltageMin"),parameterStorage->parameters.at("controlVoltageMax"));
                }
            }
            genomeSet[k].second=0;
        }

        //genome 20-25, rnd
        for(int k=20; k < 25; k++){
            for(int i=0; i < controlElectrodes; i++){
                genomeSet[k].first[i]=enhance::random_double(parameterStorage->parameters.at("controlVoltageMin"),parameterStorage->parameters.at("controlVoltageMax"));
            }
            genomeSet[k].second=0;
        }


        //mutate
        bool mutatedThisGenome;
        for(int k=0; k < 25; k++){
            mutatedThisGenome=false;
            for(int i=0; i < controlElectrodes; i++){
                if (enhance::random_double(0,1)>0.9){
                    genomeSet[k].first[i]=enhance::random_triangle(parameterStorage->parameters.at("controlVoltageMin"),genomeSet[k].first[i],parameterStorage->parameters.at("controlVoltageMax"));
                    mutatedThisGenome=true;
                }
            }
            if (mutatedThisGenome){
                genomeSet[k].second=0;
            }
        }
        // genomeSet[k].second=0 for all changed genomeSets

        // ------------ run generation ------
        std::cout<<"------------------------------ run geneartion "<<generation<< " ------------------------------"<<std::endl;
        for(int k=0; k < 25; k++){
            std::cout<<"genome: "<<k<<" voltages:";
            for(int i=0; i < controlElectrodes; i++){
                parameterStorage->electrodes[controlElectrodeIndices[i]].voltage=genomeSet[k].first[i];
                std::cout<<" "<<controlElectrodeIndices[i]<<": "<<parameterStorage->electrodes[controlElectrodeIndices[i]].voltage;
            }
            std::cout<<std::endl;

            auto startTime = chrono::steady_clock::now();

            runVoltageSetup();
            calcOptimizationEnergy();
            saveResults();
            dataFile->addData("generation",& generation);
            
            genomeSet[k].second=optEnergy;


            std::cout<<"optEnergy: "<<optEnergy    <<" fitness: ("<<fitness    <<" +- "<<fitnessUncert    <<") normedDiff: "<<normedDiff<<std::endl;
            auto endTime = chrono::steady_clock::now();
            std::cout << "time per VoltageSetup = " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count()/1000.0 << " s" << std::endl;
        }

        std::cout<<"generation "<<generation<< " done! sorted results: "<<std::endl;
        std::sort(genomeSet.begin(),genomeSet.end(),genomeComparator);
        for(int k=0; k < 25; k++){
            std::cout<<"genome: "<<k+1<<" optEnergy: "<<genomeSet[k].second<<" voltages: ";
            for(int i=0; i < controlElectrodes; i++){
                std::cout<<" "<<controlElectrodeIndices[i]<<": "<<genomeSet[k].first[i];
            }
            std::cout<<std::endl;
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