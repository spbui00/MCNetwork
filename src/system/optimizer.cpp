#include "optimizer.h"
#include "debug.h"


Optimizer::Optimizer(std::shared_ptr<ParameterStorage> parameterStorage) : 
                                                                           parameterStorage(parameterStorage),
                                                                           jobManager(parameterStorage)
{
    DEBUG_FUNC_START

    electrodeNumber=int(parameterStorage->electrodes.size());


    voltageScanPoints = parameterStorage->parameters.at("voltageScanPoints");

    for(int i=0;i<electrodeNumber;i++){
        if((i !=parameterStorage->parameters.at("outputElectrode")) & (i !=parameterStorage->parameters.at("inputElectrode1")) &(i !=parameterStorage->parameters.at("inputElectrode2"))){
            controlElectrodeIndices.push_back(i);
        }
    }
    controlElectrodeNumber = controlElectrodeIndices.size();

    voltageSets         .push_back(std::vector<double>(electrodeNumber));
    outputCurrents      .push_back(std::vector<double>(voltageScanPoints*voltageScanPoints));
    outputCurrentUncerts.push_back(std::vector<double>(voltageScanPoints*voltageScanPoints));
    voltageSets         .push_back(std::vector<double>(electrodeNumber));

    DEBUG_FUNC_END
}


/*!
  - create datafile
  - if continue mode: search for last energy point
  - start actual optimization routines
  \param startMode 0 = use voltages defined in input file, 1 = search for random start using searchForRandomStart(), 2 = continue
 */
void Optimizer::run(std::string optimizationMode, int startMode){
    DEBUG_FUNC_START

    this->optimizationMode = optimizationMode;

    //needed for genetic mode
    std::vector<std::pair<std::vector<double>,double>> startGenome = {};

    //setup data file
    std::string dataFileName = parameterStorage->workingDirecotry+ "data.hdf5";
    if (optimizationMode == "continue"){
        startMode = 2;
        //create DataFile by copying existing file
        dataFile = std::make_shared<DataFile>(dataFileName, false);

        if      (dataFile->checkDataSetExists("accepted")){
            optimizationMode = "MC";
            std::cout<< "continue MC optimization"<<std::endl;

            // search for last accepted point in datafile
            std::vector<double> accepted = * dataFile->readFullDataset("accepted"); 
            int lastAccepted=0;
            for (size_t i = accepted.size()-1; i >=0 ; i--){
                if (accepted[i] == 1){
                    lastAccepted=i;
                    break;
                }
            }
            
            //set last accepted point
            voltageSets[0]          = * dataFile->readDatasetSlice("voltages",lastAccepted); 
            outputCurrents      [0] = * dataFile->readDatasetSlice("outputCurrent",lastAccepted);
            outputCurrentUncerts[0] = * dataFile->readDatasetSlice("outputCurrentUncert",lastAccepted);

            calcOptimizationEnergy();
            
            std::cout<<"last accepted point: optEnergy: "<<optEnergy    <<" fitness: ("<<fitness    <<" +- "<<fitnessUncert    <<") normedDiff: "<<normedDiff<<std::endl;
        
            std::cout<<"voltages: "<<std::endl;
            for(int i=0;i<electrodeNumber;i++){
                if((i !=parameterStorage->parameters.at("outputElectrode")) & (i !=parameterStorage->parameters.at("inputElectrode1")) &(i !=parameterStorage->parameters.at("inputElectrode2"))){
                    std::cout<<i<<" "<<voltageSets[0][i]<<std::endl;
                }
            }




        }
        else if (dataFile->checkDataSetExists("generation")){
            optimizationMode = "genetic";

            throw std::logic_error("continuing genetic optimization - not implemented yet");

        }
        else if (dataFile->checkDataSetExists("basinNumber")){
            optimizationMode = "basinHop";

            throw std::logic_error("continuing basinHop optimization - not implemented yet");

        }
        else{
            throw std::logic_error("no dataset found in datafile that matches a started optimization -> can not continue");
        }

    }
    else{
        //creating new file
        dataFile = std::make_shared<DataFile>(dataFileName, true);
        
        dataFile->createDataset("outputCurrent",       {voltageScanPoints,voltageScanPoints});
        dataFile->createDataset("outputCurrentUncert", {voltageScanPoints,voltageScanPoints});
        dataFile->createDataset("voltages",            {electrodeNumber});
        dataFile->createDataset("fitness",             {1});
        dataFile->createDataset("fitnessUncert",       {1});
        dataFile->createDataset("optEnergy",           {1});

        //add mode specific datasets
        if      (optimizationMode == "MC"      ){ dataFile->createDataset("accepted"   ,{1}); }
        else if (optimizationMode == "genetic" ){ dataFile->createDataset("generation" ,{1}); }
        else if (optimizationMode == "basinHop"){ dataFile->createDataset("basinNumber",{1}); }
    }


    // start run
    if      (optimizationMode == "singleRun"){ singleRun();            }
    else if (optimizationMode == "MC"       ){ optimizeMC(startMode);  }
    else if (optimizationMode == "genetic"  ){ optimizeGenetic(startGenome);}
    else if (optimizationMode == "basinHop" ){ optimizeBasinHopping(); }


    DEBUG_FUNC_END
}


void Optimizer::saveResults(){
    DEBUG_FUNC_START

    dataFile->addData("outputCurrent"      ,outputCurrents      [0].data());
    dataFile->addData("outputCurrentUncert",outputCurrentUncerts[0].data());
    dataFile->addData("voltages"           ,voltageSets         [0].data());
    dataFile->addData("fitness"            ,& fitness);
    dataFile->addData("fitnessUncert"      ,& fitnessUncert);
    dataFile->addData("optEnergy"          ,& optEnergy);

    DEBUG_FUNC_END
}


void Optimizer::calcOptimizationEnergy(){
    DEBUG_FUNC_START

    int maxIndex=0,minIndex=0;

    for(int i=0; i < voltageScanPoints; i++){            
        for(int j=0; j < voltageScanPoints; j++){
            if (outputCurrents[0][i*voltageScanPoints+j]<outputCurrents[0][minIndex]){minIndex=i*voltageScanPoints+j;}
            if (outputCurrents[0][i*voltageScanPoints+j]>outputCurrents[0][maxIndex]){maxIndex=i*voltageScanPoints+j;}
        }
    }

    fitness=0;
    fitnessUncert=0;
    double normed,desiredVal,normedUncert;
    for(int i=0; i < voltageScanPoints; i++){            
        for(int j=0; j < voltageScanPoints; j++){
            normed=(outputCurrents[0][i*voltageScanPoints+j]-outputCurrents[0][minIndex])/(outputCurrents[0][maxIndex]-outputCurrents[0][minIndex]);
            normedUncert=std::sqrt(std::pow(outputCurrentUncerts[0][i*voltageScanPoints+j]/(outputCurrents[0][maxIndex]-outputCurrents[0][minIndex]),2)
                                  +std::pow((outputCurrents[0][i*voltageScanPoints+j]-outputCurrents[0][minIndex])/std::pow(outputCurrents[0][maxIndex]-outputCurrents[0][minIndex],2)*outputCurrentUncerts[0][maxIndex],2)
                                  +std::pow(((outputCurrents[0][i*voltageScanPoints+j]-outputCurrents[0][minIndex])/std::pow(outputCurrents[0][maxIndex]-outputCurrents[0][minIndex],2)-1/(outputCurrents[0][maxIndex]-outputCurrents[0][minIndex]))*outputCurrentUncerts[0][minIndex],2));
            desiredVal = desiredLogicFunction(parameterStorage->inputVoltages[i],parameterStorage->inputVoltages[j],parameterStorage->gate);
            fitness+=std::abs(normed-desiredVal);
            fitnessUncert+=normedUncert*normedUncert;
        }
    }
    fitness/=voltageScanPoints*voltageScanPoints;
    fitness=1-fitness;
    fitnessUncert=std::sqrt(fitnessUncert);
    fitnessUncert/=voltageScanPoints*voltageScanPoints;
    normedDiff= (outputCurrents[0][maxIndex]-outputCurrents[0][minIndex])/(2*std::max(std::abs(outputCurrents[0][maxIndex]),std::abs(outputCurrents[0][minIndex])));

    optEnergy = fitness - fitnessUncert*parameterStorage->parameters.at("fitnessUncertWeight") + normedDiff*parameterStorage->parameters.at("diffWeight");

    if (std::isnan(optEnergy)){
        std::cerr<<"-------------------> invalid optEnergy <-------------------"<<std::endl;
        std::cout<<"-------------------> invalid optEnergy <-------------------"<<std::endl;
        optEnergy=-INFINITY;
    }

    DEBUG_FUNC_END
}




/*!
  not performing any optimization, just runing control voltages defined in input file
 */
void Optimizer::singleRun(){
    DEBUG_FUNC_START
    std::cout<<"running fixed setup"<<std::endl;

    auto startTime = std::chrono::steady_clock::now();


    for(int i=0; i < electrodeNumber; i++){
        voltageSets[0][i] = parameterStorage->electrodes[i].voltage;
    }


    std::pair<std::vector<double>,std::vector<double>> result = jobManager.runControlVoltagesSetup(voltageSets[0]);
    outputCurrents      [0] = result.first;
    outputCurrentUncerts[0] = result.second;
    
    calcOptimizationEnergy();
    saveResults();

    std::cout<<"optEnergy: "<<optEnergy    <<" fitness: ("<<fitness    <<" +- "<<fitnessUncert    <<") normedDiff: "<<normedDiff<<std::endl;;
    auto endTime = std::chrono::steady_clock::now();
    std::cout << "time elapsed = " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count()/1000.0 << " s" << std::endl;
    
    DEBUG_FUNC_END
}


/*!
  optimize cotrol voltages using simple Monte Carlo algorithm
  \param rndStart if True: searchForRandomStart() is called to find best start point, else: voltages given in input file are used
 */
void Optimizer::optimizeMC(int startMode /*= 0*/){
    DEBUG_FUNC_START

    double accepted = 1;
    std::cout<<"running optimization - simple MC"<<std::endl;


    
    // init voltages
    if      (startMode == 0){
        for(int i=0; i < electrodeNumber; i++){
            voltageSets[0][i] = parameterStorage->electrodes[i].voltage;
        }

        std::pair<std::vector<double>,std::vector<double>> result = jobManager.runControlVoltagesSetup(voltageSets[0]);
        outputCurrents      [0] = result.first;
        outputCurrentUncerts[0] = result.second;
        calcOptimizationEnergy();
        saveResults();
        dataFile->addData("accepted",& accepted);

        std::cout<<"optEnergy: "<<optEnergy    <<" fitness: ("<<fitness    <<" +- "<<fitnessUncert    <<") normedDiff: "<<normedDiff<<std::endl;
    }
    else if (startMode == 1) {
        searchForRandomStart();
    }
    else if (startMode == 2) {
    }

    // constructive run
    double lastFitness       = fitness;
    double lastFitnessUncert = fitnessUncert;
    double lastOptEnergy     = optEnergy;
    double lastNormedDiff    = normedDiff;
    std::vector<double> lastVoltages;
    lastVoltages = voltageSets[0];

    int increaseNumber=0;
    while (optEnergy < parameterStorage->parameters.at("convergenceEnergy")){
        auto startTime = std::chrono::steady_clock::now();

        //get new random voltages
        std::cout<<"new random voltages: "<<std::endl;
        for(int i=0;i<electrodeNumber;i++){
            if((i !=parameterStorage->parameters.at("outputElectrode")) & (i !=parameterStorage->parameters.at("inputElectrode1")) &(i !=parameterStorage->parameters.at("inputElectrode2"))){
                voltageSets[0][i]=enhance::random_double(std::max(parameterStorage->parameters.at("controlVoltageMin"),voltageSets[0][i]-parameterStorage->parameters.at("maxDeltaV")),std::min(parameterStorage->parameters.at("controlVoltageMax"),voltageSets[0][i]+parameterStorage->parameters.at("maxDeltaV")));
                std::cout<<i<<" "<<voltageSets[0][i]<<std::endl;
            }
        }

        std::pair<std::vector<double>,std::vector<double>> result = jobManager.runControlVoltagesSetup(voltageSets[0]);
        outputCurrents      [0] = result.first;
        outputCurrentUncerts[0] = result.second;
        calcOptimizationEnergy();
        saveResults();

        std::cout<<"now: optEnergy: "<<optEnergy    <<" fitness: ("<<fitness    <<" +- "<<fitnessUncert    <<") normedDiff: "<<normedDiff<<
                "\nlast: optEnergy: "<<lastOptEnergy<<" fitness: ("<<lastFitness<<" +- "<<lastFitnessUncert<<") normedDiff: "<<lastNormedDiff<<std::endl;
        if((optEnergy < lastOptEnergy) & (enhance::fastExp((optEnergy-lastOptEnergy)/parameterStorage->parameters.at("MCTemp"))<enhance::random_double(0,1))){
            std::cout<<"-- not accepted --"<<std::endl;
            accepted=0;
            //swap back
            voltageSets[0] = lastVoltages;
        }
        else{
            std::cout<<"-- accepted --"<<std::endl;
            accepted=1;
            //setup for next iteration
            lastVoltages = voltageSets[0];

            lastFitness       = fitness;
            lastFitnessUncert = fitnessUncert;
            lastOptEnergy     = optEnergy;
            lastNormedDiff    = normedDiff;

        }
        dataFile->addData("accepted",& accepted);


        auto endTime = std::chrono::steady_clock::now();
        std::cout << "time per VoltageSetup = " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count()/1000.0 << " s" << std::endl;
    

        if((increaseNumber < parameterStorage->parameters.at("maxStepIncreases")) and (fitness+fitnessUncert*2)>1){
            increaseNumber++;
            parameterStorage->parameters["calcCurrentSteps"]*=2;
            std::cout<<"############ steps increased!! now: "<<parameterStorage->parameters["calcCurrentSteps"]<<" #############"<<std::endl;
        }
    }
    std::cout<<"-------------------> convergence reached <-------------------"<<std::endl;
    
    DEBUG_FUNC_END
}


void Optimizer::optimizeGenetic(std::vector<std::pair<std::vector<double>,double>> const & startGenome /* = {} */){
    DEBUG_FUNC_START

    std::cout<<"running optimization - genetic"<<std::endl;

    // lambda function needed later to sort genome set
    auto genomeComparator = []( const std::pair<std::vector<double>,double>& l, const std::pair<std::vector<double>,double>& r) { return l.second > r.second; };





    double bestFitness       = 0;
    double bestFitnessUncert = 0;


    //setup genome
    std::vector<std::pair<std::vector<double>,double>> genomeSet; //pair of control voltages and optEnergy
    for(int i=0; i < 25; i++){
        genomeSet.push_back(std::pair<std::vector<double>,double>(std::vector<double>(controlElectrodeNumber),0));
        for(int j=0; j < controlElectrodeNumber; j++){
            genomeSet[i].first[j]=enhance::random_double(parameterStorage->parameters.at("controlVoltageMin"),parameterStorage->parameters.at("controlVoltageMax"));
        }

    }

    //run first generation
    double generation=1;
    std::cout<<"------------------------------ run geneartion "<<generation<< " ------------------------------"<<std::endl;
    for(int k=0; k < 25; k++){
        std::cout<<"genome: "<<k<<" voltages:";
        for(int i=0; i < controlElectrodeNumber; i++){
            voltageSets[0][controlElectrodeIndices[i]]=genomeSet[k].first[i];
            std::cout<<" "<<controlElectrodeIndices[i]<<": "<<voltageSets[0][controlElectrodeIndices[i]];
        }
        std::cout<<std::endl;

        auto startTime = std::chrono::steady_clock::now();

        std::pair<std::vector<double>,std::vector<double>> result = jobManager.runControlVoltagesSetup(voltageSets[0]);
        outputCurrents      [0] = result.first;
        outputCurrentUncerts[0] = result.second;

        calcOptimizationEnergy();
        saveResults();
        dataFile->addData("generation",& generation);
        
        genomeSet[k].second=optEnergy;


        std::cout<<"optEnergy: "<<optEnergy    <<" fitness: ("<<fitness    <<" +- "<<fitnessUncert    <<") normedDiff: "<<normedDiff<<std::endl;
        auto endTime = std::chrono::steady_clock::now();
        std::cout << "time per VoltageSetup = " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count()/1000.0 << " s" << std::endl;
    }

    std::cout<<"generation "<<generation<< " done! sorted results: "<<std::endl;
    std::sort(genomeSet.begin(),genomeSet.end(),genomeComparator);
    for(int k=0; k < 25; k++){
        std::cout<<"genome: "<<k+1<<" optEnergy: "<<genomeSet[k].second<<" voltages: ";
        for(int i=0; i < controlElectrodeNumber; i++){
            std::cout<<" "<<controlElectrodeIndices[i]<<": "<<genomeSet[k].first[i];
        }
        std::cout<<std::endl;
    }

    int increaseNumber = 0;
    while(true){
        // ---------- setup next generation -------------
        generation++;
        bestFitness       = 0;
        bestFitnessUncert = 0;

        //first 5 genomes dont need to be changed

        //genome 6-10, rnd Bias
        for(int k=5; k < 10; k++){
            for(int i=0; i < controlElectrodeNumber; i++){
                genomeSet[k].first[i]=enhance::random_double(std::max(parameterStorage->parameters.at("controlVoltageMin"),genomeSet[k-5].first[i]-parameterStorage->parameters.at("maxDeltaV")),std::min(parameterStorage->parameters.at("controlVoltageMax"),genomeSet[k-5].first[i]+parameterStorage->parameters.at("maxDeltaV")));
            }
            genomeSet[k].second=0;
        }

        //genome 11-15, crossover
        for(int k=10; k < 15; k++){
            for(int i=0; i < controlElectrodeNumber; i++){
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
            for(int i=0; i < controlElectrodeNumber; i++){
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
            for(int i=0; i < controlElectrodeNumber; i++){
                genomeSet[k].first[i]=enhance::random_double(parameterStorage->parameters.at("controlVoltageMin"),parameterStorage->parameters.at("controlVoltageMax"));
            }
            genomeSet[k].second=0;
        }


        //mutate
        bool mutatedThisGenome;
        for(int k=0; k < 25; k++){
            mutatedThisGenome=false;
            for(int i=0; i < controlElectrodeNumber; i++){
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
            for(int i=0; i < controlElectrodeNumber; i++){
                voltageSets[0][controlElectrodeIndices[i]]=genomeSet[k].first[i];
                std::cout<<" "<<controlElectrodeIndices[i]<<": "<<voltageSets[0][controlElectrodeIndices[i]];
            }
            std::cout<<std::endl;

            auto startTime = std::chrono::steady_clock::now();

            std::pair<std::vector<double>,std::vector<double>> result = jobManager.runControlVoltagesSetup(voltageSets[0]);
            outputCurrents      [0] = result.first;
            outputCurrentUncerts[0] = result.second;

            calcOptimizationEnergy();
            saveResults();
            dataFile->addData("generation",& generation);
            
            genomeSet[k].second=optEnergy;


            std::cout<<"optEnergy: "<<optEnergy    <<" fitness: ("<<fitness    <<" +- "<<fitnessUncert    <<") normedDiff: "<<normedDiff<<std::endl;
            auto endTime = std::chrono::steady_clock::now();
            std::cout << "time per VoltageSetup = " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count()/1000.0 << " s" << std::endl;

            if (fitness > bestFitness){
                bestFitness = fitness;
                bestFitnessUncert = fitnessUncert;
            }
        }

        std::cout<<"generation "<<generation<< " done! sorted results: "<<std::endl;
        std::sort(genomeSet.begin(),genomeSet.end(),genomeComparator);
        for(int k=0; k < 25; k++){
            std::cout<<"genome: "<<k+1<<" optEnergy: "<<genomeSet[k].second<<" voltages: ";
            for(int i=0; i < controlElectrodeNumber; i++){
                std::cout<<" "<<controlElectrodeIndices[i]<<": "<<genomeSet[k].first[i];
            }
            std::cout<<std::endl;
        }


        if((increaseNumber < parameterStorage->parameters.at("maxStepIncreases")) and (bestFitness+bestFitnessUncert*2)>1){
            increaseNumber++;
            parameterStorage->parameters["calcCurrentSteps"]*=2;
            std::cout<<"############ steps increased!! now: "<<parameterStorage->parameters.at("calcCurrentSteps")<<" #############"<<std::endl;
        }
        if (genomeSet[0].second > parameterStorage->parameters.at("convergenceEnergy")){
            break;
        }

    }    

    DEBUG_FUNC_END
}

void Optimizer::optimizeBasinHopping(bool rndStart /*= false*/){
    DEBUG_FUNC_START

    throw std::logic_error("basin hopping - not implemented yet");

    std::cout<<"running optimization - basin hopping"<<std::endl;


    


    // init random voltages
    if (rndStart){
        searchForRandomStart();
    }
    else{
        for(int i=0; i < electrodeNumber; i++){
            voltageSets[0][i] = parameterStorage->electrodes[i].voltage;
        }
    }

    // constructive run
    double accepted          = -1;
    double lastFitness       = fitness;
    double lastFitnessUncert = fitnessUncert;
    double lastOptEnergy     = optEnergy;
    double lastNormedDiff    = normedDiff;
    std::vector<double> lastVoltages = voltageSets[0];

    for (size_t i = 0; i < 1000000; i++){
        auto startTime = std::chrono::steady_clock::now();

        //get new random voltages
        std::cout<<"new random voltages: "<<std::endl;
        for(int i=0;i<electrodeNumber;i++){
            if((i !=parameterStorage->parameters.at("outputElectrode")) & (i !=parameterStorage->parameters.at("inputElectrode1")) &(i !=parameterStorage->parameters.at("inputElectrode2"))){
                voltageSets[0][i]=enhance::random_double(std::max(parameterStorage->parameters.at("controlVoltageMin"),voltageSets[0][i]-parameterStorage->parameters.at("basinDeltaV")),std::min(parameterStorage->parameters.at("controlVoltageMax"),voltageSets[0][i]+parameterStorage->parameters.at("basinDeltaV")));
                std::cout<<i<<" "<<voltageSets[0][i]<<std::endl;
            }
        }


        optimizeGradient(i);

        std::cout<<"now: optEnergy: "<<optEnergy    <<" fitness: ("<<fitness    <<" +- "<<fitnessUncert    <<") normedDiff: "<<normedDiff<<
                "\nlast: optEnergy: "<<lastOptEnergy<<" fitness: ("<<lastFitness<<" +- "<<lastFitnessUncert<<") normedDiff: "<<lastNormedDiff<<std::endl;
        if((optEnergy < lastOptEnergy) & (enhance::fastExp((optEnergy-lastOptEnergy)/parameterStorage->parameters.at("basinTemp"))<enhance::random_double(0,1))){
            std::cout<<"-- not accepted --"<<std::endl;
            accepted=0;
            //swap back
            voltageSets[0] = lastVoltages;
        }
        else{
            std::cout<<"-- accepted --"<<std::endl;
            accepted=1;
            //setup for next iteration
            lastVoltages = voltageSets[0];

            lastFitness       = fitness;
            lastFitnessUncert = fitnessUncert;
            lastOptEnergy     = optEnergy;
            lastNormedDiff    = normedDiff;

        }

        auto endTime = std::chrono::steady_clock::now();
        std::cout << "time per VoltageSetup = " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count()/1000.0 << " s" << std::endl;
    

    }
    

    DEBUG_FUNC_END
}

void Optimizer::optimizeGradient(int basinNumber){
    DEBUG_FUNC_START


    std::vector<double> gradient(controlElectrodeNumber);

    std::vector<double> lastVoltages = voltageSets[0];
    double lastOptEnergy = optEnergy;
    int gradientComponentSign = 1; //if V+deltaV to cals gradient exceeds controlVoltageMax,  V-deltaV ist used instead. info ist stored in this variable
    double stepWidth = 1;
    while (true){
        //calc Gradient
        for (size_t i = 0; i < controlElectrodeNumber; i++){
            voltageSets[0] = lastVoltages;
            if (voltageSets[0][controlElectrodeIndices[i]] + parameterStorage->parameters.at("gradDeltaV") <  parameterStorage->parameters.at("controlVoltageMax")){
                gradientComponentSign = 1;
                voltageSets[0][controlElectrodeIndices[i]] +=  parameterStorage->parameters.at("gradDeltaV");
            }
            else{
                gradientComponentSign = -1;
                voltageSets[0][controlElectrodeIndices[i]] -=  parameterStorage->parameters.at("gradDeltaV");
            }

            std::pair<std::vector<double>,std::vector<double>> result = jobManager.runControlVoltagesSetup(voltageSets[0]);
            outputCurrents      [0] = result.first;
            outputCurrentUncerts[0] = result.second;

            calcOptimizationEnergy();

            gradient[i] = gradientComponentSign*(lastOptEnergy-optEnergy)/parameterStorage->parameters.at("gradDeltaV");
        }
        
        //move Grad Step


        if (false){ //converged
            break;
        }
    }
    

    DEBUG_FUNC_END
}

/*!
    generates "rndStartPoints" random control voltage points and sets voltageSets[0] and optEnergy to the best result
*/
void Optimizer::searchForRandomStart(){
    DEBUG_FUNC_START

    std::cout<<"------ searching for start point ------"<<std::endl;

    std::vector<int> controlElectrodeIndices;
    for(int i=0;i<electrodeNumber;i++){
        if((i !=parameterStorage->parameters.at("outputElectrode")) & (i !=parameterStorage->parameters.at("inputElectrode1")) &(i !=parameterStorage->parameters.at("inputElectrode2"))){
            controlElectrodeIndices.push_back(i);
        }
    }
    int controlElectrodeNumber=controlElectrodeIndices.size();
    
    std::vector<std::pair<std::vector<double>,double>> rndStartPoints; //pair of control voltages and optEnergy
    auto rndStartPointsComparator = []( const std::pair<std::vector<double>,double>& l, const std::pair<std::vector<double>,double>& r) { return l.second > r.second; }; //needed to find best start point

    for(int k=0; k < parameterStorage->parameters.at("rndStartPoints"); k++){
        rndStartPoints.push_back(std::pair<std::vector<double>,double>(std::vector<double>(controlElectrodeNumber),0));
        for(int j=0; j < controlElectrodeNumber; j++){
            rndStartPoints[k].first[j]=enhance::random_double(parameterStorage->parameters.at("controlVoltageMin"),parameterStorage->parameters.at("controlVoltageMax"));
        }
    }

    //run rnd start candidates
    for(int k=0; k < parameterStorage->parameters.at("rndStartPoints"); k++){
        std::cout<<"rnd start point: "<<k<<" voltages:";
        for(int i=0; i < controlElectrodeNumber; i++){
            voltageSets[0][controlElectrodeIndices[i]]=rndStartPoints[k].first[i];
            std::cout<<" "<<controlElectrodeIndices[i]<<": "<<voltageSets[0][controlElectrodeIndices[i]];
        }
        std::cout<<std::endl;

        auto startTime = std::chrono::steady_clock::now();

        std::pair<std::vector<double>,std::vector<double>> result = jobManager.runControlVoltagesSetup(voltageSets[0]);
        outputCurrents      [0] = result.first;
        outputCurrentUncerts[0] = result.second;

        calcOptimizationEnergy();
        saveResults();
        if (optimizationMode == "MC"){
            int a = 1;
            dataFile->addData("accepted",& a);
        }
        else if (optimizationMode == "basinHop"){
            int a = -1;
            dataFile->addData("basinNumber",& a);
        }
        
        rndStartPoints[k].second=optEnergy;


        std::cout<<"optEnergy: "<<optEnergy    <<" fitness: ("<<fitness    <<" +- "<<fitnessUncert    <<") normedDiff: "<<normedDiff<<std::endl;
        auto endTime = std::chrono::steady_clock::now();
        std::cout << "time per VoltageSetup = " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count()/1000.0 << " s" << std::endl;
    }

    std::cout<<"start search done! sorted results: "<<std::endl;
    std::sort(rndStartPoints.begin(),rndStartPoints.end(),rndStartPointsComparator);
    for(int k=0; k < parameterStorage->parameters.at("rndStartPoints"); k++){
        std::cout<<"startPoint: "<<k+1<<" optEnergy: "<<rndStartPoints[k].second<<" voltages: ";
        for(int i=0; i < controlElectrodeNumber; i++){
            std::cout<<" "<<controlElectrodeIndices[i]<<": "<<rndStartPoints[k].first[i];
        }
        std::cout<<std::endl;
    }

    //set best start
    for(int i=0; i < controlElectrodeNumber; i++){
        voltageSets[0][controlElectrodeIndices[i]]=rndStartPoints[0].first[i];
    }
    optEnergy = rndStartPoints[0].second;
    

    DEBUG_FUNC_END
}


bool Optimizer::desiredLogicFunction(double val1, double val2, std::string gate){
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
        throw std::runtime_error("logic gate not found");
    }

    DEBUG_FUNC_END
}
