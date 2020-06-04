#include "optimizer.h"
#include "debug.h"

Optimizer::Optimizer(std::shared_ptr<ParameterStorage> parameterStorage)
    : parameterStorage(parameterStorage)
    , jobManager(parameterStorage)
{
    DEBUG_FUNC_START

    electrodeNumber = int(parameterStorage->electrodes.size());

    voltageScanPoints = parameterStorage->parameters.at("voltageScanPoints");

    for (int i = 0; i < electrodeNumber; i++) {
        if ((i != parameterStorage->parameters.at("outputElectrode")) & (i != parameterStorage->parameters.at("inputElectrode1")) & (i != parameterStorage->parameters.at("inputElectrode2"))) {
            controlElectrodeIndices.push_back(i);
        }
    }
    controlElectrodeNumber = controlElectrodeIndices.size();

    voltageEnergySets.push_back(std::pair<std::vector<double>, double>(
        std::vector<double>(electrodeNumber), 0));

    DEBUG_FUNC_END
}

/*!
  - create datafile
  - start actual optimization routines
  \param startMode 0 = use voltages defined in input file, 1 = search for random
  start using searchForRandomStart(), 2 = continue
 */
void Optimizer::run(std::string optimizationMode, int startMode)
{
    DEBUG_FUNC_START

    this->optimizationMode = optimizationMode;

    if (optimizationMode == "continue") {
        startMode = 2;
        continueSimulation();
    } else {
        // creating new file
        dataFile = std::make_shared<DataFile>(
            parameterStorage->workingDirecotry + "data.hdf5", true);

        dataFile->createDataset("outputCurrent",
            { voltageScanPoints, voltageScanPoints });
        dataFile->createDataset("outputCurrentUncert",
            { voltageScanPoints, voltageScanPoints });
        dataFile->createDataset("voltages", { electrodeNumber });
        dataFile->createDataset("fitness", { 1 });
        dataFile->createDataset("fitnessUncert", { 1 });
        dataFile->createDataset("optEnergy", { 1 });

        // add mode specific datasets
        if (this->optimizationMode == "MC") {
            dataFile->createDataset("accepted", { 1 });
        } // -1 = start, 0 = false, 1 = true
        else if (this->optimizationMode == "genetic") {
            dataFile->createDataset("generation", { 1 });
        } // generation number
        else if (this->optimizationMode == "basinHop") {
            dataFile->createDataset("basinAccepted", { 1 });
        } // -1 = start, 0 = flase, 1 =true, 2 = basin jump
    }

    // start run
    if (this->optimizationMode == "singleRun") {
        singleRun(startMode);
    } else if (this->optimizationMode == "MC") {
        optimizeMC(startMode);
    } else if (this->optimizationMode == "genetic") {
        optimizeGenetic(startMode);
    } else if (this->optimizationMode == "basinHop") {
        optimizeBasinHopping(startMode);
    }

    DEBUG_FUNC_END
}

/*!
    is called by Optimizer::run in case of "--continue" option set. setup
   optimization parameters for continuation.
 */
void Optimizer::continueSimulation()
{
    DEBUG_FUNC_START

    // create DataFile by copying existing file
    dataFile = std::make_shared<DataFile>(
        parameterStorage->workingDirecotry + "data.hdf5", false);

    // check existing datasets in datafile to determine optimization mode
    if (dataFile->checkDataSetExists("accepted")) {
        optimizationMode = "MC";
        std::cout << "continue MC optimization" << std::endl;

        // search for last accepted point in datafile
        std::vector<double> accepted = *dataFile->readFullDataset("accepted");
        iteration = accepted.size();

        int lastAccepted = 0;
        for (int i = accepted.size() - 1; i >= 0; i--) {
            if (accepted[i] == 1) {
                lastAccepted = i;
                break;
            }
        }

        // set last accepted point
        voltageEnergySets[0].first = *dataFile->readDatasetSlice("voltages", lastAccepted);
        outputCurrents = *dataFile->readDatasetSlice("outputCurrent", lastAccepted);
        outputCurrentUncerts = *dataFile->readDatasetSlice("outputCurrentUncert", lastAccepted);

        calcOptimizationEnergy();
        voltageEnergySets[0].second = optEnergy;

        std::cout << "last accepted point: optEnergy: "
                  << voltageEnergySets[0].second << " fitness: (" << fitness
                  << " +- " << fitnessUncert << ") normedDiff: " << normedDiff
                  << std::endl;

        std::cout << "voltages: " << std::endl;
        for (int i = 0; i < controlElectrodeNumber; i++) {
            std::cout << controlElectrodeIndices[i] << " "
                      << voltageEnergySets[0].first[controlElectrodeIndices[i]]
                      << std::endl;
        }

    } else if (dataFile->checkDataSetExists("generation")) {
        optimizationMode = "genetic";
        std::cout << "continue genetic optimization" << std::endl;

        while (voltageEnergySets.size() < 25) {
            voltageEnergySets.push_back(std::pair<std::vector<double>, double>(
                std::vector<double>(electrodeNumber), 0));
        }

        std::vector<double> generations = *dataFile->readFullDataset("generation");
        double lastGeneration = generations.back();
        size_t generationSize = 0, lastGenerationStartIdx = 0;

        for (int i = generations.size() - 1; i >= 0; i--) {
            if (lastGeneration == generations[i]) {
                generationSize++;
                if (generationSize == 25) {
                    lastGenerationStartIdx = i;
                    iteration = lastGeneration; // iteration used as buffer for generation here
                    std::cout << "last full generation found! gen: " << lastGeneration
                              << " at index " << lastGenerationStartIdx << std::endl;
                    break;
                }
            } else {
                // deleting incomplete generation
                dataFile->shrinkDataset("fitness", generationSize);
                dataFile->shrinkDataset("fitnessUncert", generationSize);
                dataFile->shrinkDataset("generation", generationSize);
                dataFile->shrinkDataset("optEnergy", generationSize);
                dataFile->shrinkDataset("outputCurrent", generationSize);
                dataFile->shrinkDataset("outputCurrentUncert", generationSize);
                dataFile->shrinkDataset("voltages", generationSize);

                generationSize = 1;
                lastGeneration = generations[i];
            }
        }

        if (generationSize != 25) {
            throw std::logic_error("no full generation found -> can not continue");
        }

        // read top 5
        for (size_t k = 0; k < 25; k++) {
            voltageEnergySets[k].first = *dataFile->readDatasetSlice("voltages", lastGenerationStartIdx + k);
            voltageEnergySets[k].second = (*dataFile->readDatasetSlice(
                "optEnergy", lastGenerationStartIdx + k))[0];
        }

    } else if (dataFile->checkDataSetExists("basinAccepted")) {
        optimizationMode = "basinHop";
        std::cout << "continue basinHop optimization" << std::endl;

        while (voltageEnergySets.size() < 4) {
            voltageEnergySets.push_back(std::pair<std::vector<double>, double>(
                std::vector<double>(electrodeNumber), 0));
        }

        std::vector<double> accepted = *dataFile->readFullDataset("basinAccepted");
        std::vector<double> optEnergy = *dataFile->readFullDataset("optEnergy");
        int index = 0;
        iteration = accepted.size();

        // search for last accepted point
        for (size_t i = accepted.size() - 1; i >= 0; i--) {
            if (accepted[i] == 1) {
                index = i;
                break;
            }
        }
        voltageEnergySets[0].first = *dataFile->readDatasetSlice("voltages", index);
        voltageEnergySets[0].second = optEnergy[index];

        // search for current basin best
        double best = -INFINITY;
        for (size_t i = accepted.size() - 1; i >= 0; i--) {
            if (optEnergy[i] > best) {
                best = optEnergy[i];
                index = i;
            }
            if (accepted[i] == 2 | accepted[i] == 3) {
                break;
            }
        }
        voltageEnergySets[2].first = *dataFile->readDatasetSlice("voltages", index);
        voltageEnergySets[2].second = best;
        lastIterationIncrease = index;

        // search for last accepted basin best
        // first search for last accepted basin
        best = -INFINITY;
        index = 0;
        for (size_t i = accepted.size() - 1; i >= 0; i--) {
            if (accepted[i] == 3) {
                index = i;
                break;
            }
        }
        // search for best energy in basin
        for (int i = index - 1; i >= 0; i--) { // starting at last basin
            if (optEnergy[i] > best) {
                best = optEnergy[i];
                index = i;
            }
            if (accepted[i] == 2 | accepted[i] == 3) {
                break;
            }
        }
        voltageEnergySets[3].first = *dataFile->readDatasetSlice("voltages", index);
        voltageEnergySets[3].second = best;

        std::cout << "last accepted:      " << voltageEnergySets[0].second
                  << std::endl;
        std::cout << "current basin best: " << voltageEnergySets[2].second
                  << std::endl;
        std::cout << "last    basin best: " << voltageEnergySets[3].second
                  << std::endl;

    } else {
        throw std::logic_error("no dataset found in datafile that matches a "
                               "started optimization -> can not continue");
    }

    DEBUG_FUNC_END
}

/*!
    saves following datsets to dataFile: Optimizer::outputCurrent,
   Optimizer::outputCurrentUncert, Optimizer::voltages, Optimizer::fitness,
   Optimizer::fitnessUncert, oOptimizer::ptEnergy \param index determines which
   set stored in Optimizer::voltageEnergySets shall be saved. default = 0.
 */
void Optimizer::saveResults(size_t index /* = 0 */)
{
    DEBUG_FUNC_START

    dataFile->addData("outputCurrent", outputCurrents.data());
    dataFile->addData("outputCurrentUncert", outputCurrentUncerts.data());
    dataFile->addData("voltages", voltageEnergySets[index].first.data());
    dataFile->addData("fitness", &fitness);
    dataFile->addData("fitnessUncert", &fitnessUncert);
    dataFile->addData("optEnergy", &voltageEnergySets[index].second);

    DEBUG_FUNC_END
}

/*!
    calcOptimizationEnergy using values stored in Optimizer::outputCurrents and
   Optimizer::outputCurrentUncerts. saves output to Optimizer::fitness,
   Optimizer::fitnessUncert, Optimizer::normedDiff, Optimizer::optEnergy \param
   index determines which set stored in voltageEnergySets shall be saved.
   default = 0.
 */
void Optimizer::calcOptimizationEnergy()
{
    DEBUG_FUNC_START

    int maxIndex = 0, minIndex = 0;

    for (int i = 0; i < voltageScanPoints; i++) {
        for (int j = 0; j < voltageScanPoints; j++) {
            if (outputCurrents[i * voltageScanPoints + j] < outputCurrents[minIndex]) {
                minIndex = i * voltageScanPoints + j;
            }
            if (outputCurrents[i * voltageScanPoints + j] > outputCurrents[maxIndex]) {
                maxIndex = i * voltageScanPoints + j;
            }
        }
    }

    fitness = 0;
    fitnessUncert = 0;
    double normed, desiredVal, normedUncert;
    for (int i = 0; i < voltageScanPoints; i++) {
        for (int j = 0; j < voltageScanPoints; j++) {
            normed = (outputCurrents[i * voltageScanPoints + j] - outputCurrents[minIndex]) / (outputCurrents[maxIndex] - outputCurrents[minIndex]);
            normedUncert = std::sqrt(
                std::pow(outputCurrentUncerts[i * voltageScanPoints + j] / (outputCurrents[maxIndex] - outputCurrents[minIndex]),
                    2)
                + std::pow(
                    (outputCurrents[i * voltageScanPoints + j] - outputCurrents[minIndex]) / std::pow(outputCurrents[maxIndex] - outputCurrents[minIndex], 2) * outputCurrentUncerts[maxIndex],
                    2)
                + std::pow(
                    ((outputCurrents[i * voltageScanPoints + j] - outputCurrents[minIndex]) / std::pow(outputCurrents[maxIndex] - outputCurrents[minIndex], 2) - 1 / (outputCurrents[maxIndex] - outputCurrents[minIndex])) * outputCurrentUncerts[minIndex],
                    2));
            desiredVal = desiredLogicFunction(parameterStorage->inputVoltages[i],
                parameterStorage->inputVoltages[j],
                parameterStorage->gate);
            fitness += std::abs(normed - desiredVal);
            fitnessUncert += normedUncert * normedUncert;
        }
    }
    fitness /= voltageScanPoints * voltageScanPoints;
    fitness = 1 - fitness;
    fitnessUncert = std::sqrt(fitnessUncert);
    fitnessUncert /= voltageScanPoints * voltageScanPoints;
    normedDiff = (outputCurrents[maxIndex] - outputCurrents[minIndex]) / (2 * std::max(std::abs(outputCurrents[maxIndex]), std::abs(outputCurrents[minIndex])));

    optEnergy = fitness - fitnessUncert * parameterStorage->parameters.at("fitnessUncertWeight") + normedDiff * parameterStorage->parameters.at("diffWeight");

    if (std::isnan(optEnergy)) {
        std::cerr << "-------------------> invalid optEnergy <-------------------"
                  << std::endl;
        std::cout << "-------------------> invalid optEnergy <-------------------"
                  << std::endl;
        optEnergy = -INFINITY;
    }

    DEBUG_FUNC_END
}

/*!
  not performing any optimization, just runing control voltages defined in input
  file
 */
void Optimizer::singleRun(size_t startMode)
{
    DEBUG_FUNC_START
    std::cout << "running fixed setup" << std::endl;

    auto startTime = std::chrono::steady_clock::now();

    if (startMode == 0) { // voltages given in input
        for (int i = 0; i < electrodeNumber; i++) {
            voltageEnergySets[0].first[i] = parameterStorage->electrodes[i].voltage;
        }
    } else if (startMode == 1) { // rnd point
        for (int i = 0; i < electrodeNumber; i++) {
            voltageEnergySets[0].first[i] = enhance::random_double(
                parameterStorage->parameters.at("controlVoltageMin"),
                parameterStorage->parameters.at("controlVoltageMax"));
        }
    }

    std::pair<std::vector<double>, std::vector<double>> result = jobManager.runControlVoltagesSetup(voltageEnergySets[0].first);
    outputCurrents = result.first;
    outputCurrentUncerts = result.second;

    calcOptimizationEnergy();
    voltageEnergySets[0].second = optEnergy;
    saveResults(0);

    std::cout << "optEnergy: " << voltageEnergySets[0].second << " fitness: ("
              << fitness << " +- " << fitnessUncert
              << ") normedDiff: " << normedDiff << std::endl;
    ;
    auto endTime = std::chrono::steady_clock::now();
    std::cout << "time elapsed = "
              << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime)
                     .count()
            / 1000.0
              << " s" << std::endl;

    DEBUG_FUNC_END
}

/*!
  optimize cotrol voltages using simple Monte Carlo algorithm
  \param rndStart 0: searchForRandomStart() is called to find best start point
  \n 1: voltages given in input file are used \n 2: used in continue mode,
  voltages have been set before
 */
void Optimizer::optimizeMC(size_t startMode /*= 0*/)
{
    DEBUG_FUNC_START
    // voltageEnergySets positioning: 0: current point. 1: last point
    while (voltageEnergySets.size() < 2) {
        voltageEnergySets.push_back(std::pair<std::vector<double>, double>(
            std::vector<double>(electrodeNumber), 0));
    }

    double accepted = -1;
    std::cout << "running optimization - simple MC" << std::endl;

    // init voltages
    if (startMode == 0) { // voltages given in input
        for (int i = 0; i < electrodeNumber; i++) {
            voltageEnergySets[0].first[i] = parameterStorage->electrodes[i].voltage;
        }

        std::pair<std::vector<double>, std::vector<double>> result = jobManager.runControlVoltagesSetup(voltageEnergySets[0].first);
        outputCurrents = result.first;
        outputCurrentUncerts = result.second;
        calcOptimizationEnergy();
        voltageEnergySets[0].second = optEnergy;
        saveResults(0);
        dataFile->addData("accepted", &accepted);

        std::cout << "optEnergy: " << voltageEnergySets[0].second << " fitness: ("
                  << fitness << " +- " << fitnessUncert
                  << ") normedDiff: " << normedDiff << std::endl;
    } else if (startMode == 1) { // rnd start
        searchForRandomStart();
    } else if (startMode == 2) { // start set before (continue mode)
    }

    // constructive run
    double lastFitness = fitness;
    double lastFitnessUncert = fitnessUncert;
    double lastNormedDiff = normedDiff;

    voltageEnergySets[1] = voltageEnergySets[0];

    int increaseNumber = 0;
    while (optEnergy < parameterStorage->parameters.at("convergenceEnergy") & iteration < parameterStorage->parameters.at("maxIterations")) {
        iteration++;
        auto startTime = std::chrono::steady_clock::now();

        // get new random voltages
        std::cout << "new random voltages: " << std::endl;
        for (int i = 0; i < controlElectrodeNumber; i++) {
            voltageEnergySets[0].first[controlElectrodeIndices[i]] = enhance::random_double(
                std::max(parameterStorage->parameters.at("controlVoltageMin"),
                    voltageEnergySets[0].first[controlElectrodeIndices[i]] - parameterStorage->parameters.at("maxDeltaV")),
                std::min(parameterStorage->parameters.at("controlVoltageMax"),
                    voltageEnergySets[0].first[controlElectrodeIndices[i]] + parameterStorage->parameters.at("maxDeltaV")));
            std::cout << controlElectrodeIndices[i] << " "
                      << voltageEnergySets[0].first[controlElectrodeIndices[i]]
                      << std::endl;
        }

        // run
        std::pair<std::vector<double>, std::vector<double>> result = jobManager.runControlVoltagesSetup(voltageEnergySets[0].first);
        outputCurrents = result.first;
        outputCurrentUncerts = result.second;
        calcOptimizationEnergy();
        voltageEnergySets[0].second = optEnergy;
        saveResults(0);

        std::cout << "iteration " << iteration << std::endl;
        std::cout << "now:  optEnergy: " << voltageEnergySets[0].second
                  << " fitness: (" << fitness << " +- " << fitnessUncert
                  << ") normedDiff: " << normedDiff
                  << "\nlast: optEnergy: " << voltageEnergySets[1].second
                  << " fitness: (" << lastFitness << " +- " << lastFitnessUncert
                  << ") normedDiff: " << lastNormedDiff << std::endl;
        if ((voltageEnergySets[0].second < voltageEnergySets[1].second) & (enhance::fastExp((voltageEnergySets[0].second - voltageEnergySets[1].second) / parameterStorage->parameters.at("MCTemp")) < enhance::random_double(0, 1))) {
            std::cout << "-- not accepted --" << std::endl;
            accepted = 0;
            // swap back
            voltageEnergySets[0].first = voltageEnergySets[1].first;
        } else {
            std::cout << "-- accepted --" << std::endl;
            accepted = 1;
            // setup for next iteration
            voltageEnergySets[1] = voltageEnergySets[0];

            lastFitness = fitness;
            lastFitnessUncert = fitnessUncert;
            lastNormedDiff = normedDiff;
        }
        dataFile->addData("accepted", &accepted);

        auto endTime = std::chrono::steady_clock::now();
        std::cout << "time per VoltageSetup = "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(
                         endTime - startTime)
                         .count()
                / 1000.0
                  << " s" << std::endl;

        if ((increaseNumber < parameterStorage->parameters.at("maxStepIncreases")) and (fitness + fitnessUncert * 2) > 1) {
            increaseNumber++;
            parameterStorage->parameters["calcCurrentSteps"] *= 2;
            std::cout << "############ steps increased!! now: "
                      << parameterStorage->parameters["calcCurrentSteps"]
                      << " #############" << std::endl;
        }
    }
    std::cout << "-------------------> optimization stopped <-------------------"
              << std::endl;
    std::cout << "iteration: " << iteration << " optEnergy: " << optEnergy
              << std::endl;

    DEBUG_FUNC_END
}

/*!
  optimize cotrol voltages using genetic algorithm
  \param rndStart 0: searchForRandomStart() is called to find best start point
  \n 1: voltages given in input file are used \n 2: used in continue mode,
  voltages have been set before
 */
void Optimizer::optimizeGenetic(size_t startMode /* = 0 */)
{
    DEBUG_FUNC_START
    // voltageEnergySets positioning: 25 genomes
    while (voltageEnergySets.size() < 25) {
        voltageEnergySets.push_back(std::pair<std::vector<double>, double>(
            std::vector<double>(electrodeNumber), 0));
    }

    std::cout << "running optimization - genetic" << std::endl;

    // lambda function needed later to sort genome set
    auto genomeComparator = [](const std::pair<std::vector<double>, double>& l,
                                const std::pair<std::vector<double>, double>& r) {
        return l.second > r.second;
    };

    double bestFitness = 0;
    double bestFitnessUncert = 0;
    double generation = 1;

    if (startMode == 0 or startMode == 1) { // new run
        // setup genome
        for (int k = 0; k < 25; k++) {
            for (int i = 0; i < controlElectrodeNumber; i++) {
                voltageEnergySets[k].first[controlElectrodeIndices[i]] = enhance::random_double(
                    parameterStorage->parameters.at("controlVoltageMin"),
                    parameterStorage->parameters.at("controlVoltageMax"));
            }
        }

        // run first generation
        std::cout << "------------------------------ run geneartion " << generation
                  << " ------------------------------" << std::endl;
        for (int k = 0; k < 25; k++) {
            std::cout << "genome: " << k << " voltages:";
            for (int i = 0; i < controlElectrodeNumber; i++) {
                std::cout << " " << controlElectrodeIndices[i] << ": "
                          << voltageEnergySets[k].first[controlElectrodeIndices[i]];
            }
            std::cout << std::endl;

            auto startTime = std::chrono::steady_clock::now();

            std::pair<std::vector<double>, std::vector<double>> result = jobManager.runControlVoltagesSetup(voltageEnergySets[k].first);
            outputCurrents = result.first;
            outputCurrentUncerts = result.second;

            calcOptimizationEnergy();
            voltageEnergySets[k].second = optEnergy;
            saveResults(k);
            dataFile->addData("generation", &generation);

            std::cout << "optEnergy: " << voltageEnergySets[k].second << " fitness: ("
                      << fitness << " +- " << fitnessUncert
                      << ") normedDiff: " << normedDiff << std::endl;
            auto endTime = std::chrono::steady_clock::now();
            std::cout << "time per VoltageSetup = "
                      << std::chrono::duration_cast<std::chrono::milliseconds>(
                             endTime - startTime)
                             .count()
                    / 1000.0
                      << " s" << std::endl;
        }

        std::cout << "generation " << generation
                  << " done! sorted results: " << std::endl;
        std::sort(voltageEnergySets.begin(), voltageEnergySets.end(),
            genomeComparator);
        for (int k = 0; k < 25; k++) {
            std::cout << "genome: " << k + 1
                      << " optEnergy: " << voltageEnergySets[k].second
                      << " voltages: ";
            for (int i = 0; i < controlElectrodeNumber; i++) {
                std::cout << " " << controlElectrodeIndices[i] << ": "
                          << voltageEnergySets[k].first[controlElectrodeIndices[i]];
            }
            std::cout << std::endl;
        }
    } else if (startMode == 2) { // start set before (continue mode)
        generation = iteration;
        std::cout << "last generation " << generation << " was " << std::endl;
        std::sort(voltageEnergySets.begin(), voltageEnergySets.end(),
            genomeComparator);
        for (int k = 0; k < 25; k++) {
            std::cout << "genome: " << k + 1
                      << " optEnergy: " << voltageEnergySets[k].second
                      << " voltages: ";
            for (int i = 0; i < controlElectrodeNumber; i++) {
                std::cout << " " << controlElectrodeIndices[i] << ": "
                          << voltageEnergySets[k].first[controlElectrodeIndices[i]];
            }
            std::cout << std::endl;
        }
    }

    int increaseNumber = 0;
    while (true) {
        // ---------- setup next generation -------------
        generation++;
        bestFitness = 0;
        bestFitnessUncert = 0;

        // first 5 genomes dont need to be changed

        // genome 6-10, rnd Bias
        for (int k = 5; k < 10; k++) {
            for (int i = 0; i < controlElectrodeNumber; i++) {
                voltageEnergySets[k].first[controlElectrodeIndices[i]] = enhance::random_double(
                    std::max(
                        parameterStorage->parameters.at("controlVoltageMin"),
                        voltageEnergySets[k - 5].first[controlElectrodeIndices[i]] - parameterStorage->parameters.at("maxDeltaV")),
                    std::min(
                        parameterStorage->parameters.at("controlVoltageMax"),
                        voltageEnergySets[k - 5].first[controlElectrodeIndices[i]] + parameterStorage->parameters.at("maxDeltaV")));
            }
            voltageEnergySets[k].second = 0;
        }

        // genome 11-15, crossover
        for (int k = 10; k < 15; k++) {
            for (int i = 0; i < controlElectrodeNumber; i++) {
                if (enhance::random_double(0, 1) > 0.5) {
                    voltageEnergySets[k].first[controlElectrodeIndices[i]] = voltageEnergySets[k - 10].first[controlElectrodeIndices[i]];
                } else {
                    voltageEnergySets[k].first[controlElectrodeIndices[i]] = voltageEnergySets[k - 9]
                                                                                 .first[controlElectrodeIndices[i]]; ///<<<< bug here: if k=14,
                    ///<k=5 is taken, but
                    ///<voltageEnergySets[5]
                    ///<was overwritten before
                }
            }
            voltageEnergySets[k].second = 0;
        }

        // genome 16-20, rnd crossover
        for (int k = 15; k < 20; k++) {
            for (int i = 0; i < controlElectrodeNumber; i++) {
                if (enhance::random_double(0, 1) > 0.5) {
                    voltageEnergySets[k].first[controlElectrodeIndices[i]] = voltageEnergySets[k - 15].first[controlElectrodeIndices[i]];
                } else {
                    voltageEnergySets[k].first[controlElectrodeIndices[i]] = enhance::random_double(
                        parameterStorage->parameters.at("controlVoltageMin"),
                        parameterStorage->parameters.at("controlVoltageMax"));
                }
            }
            voltageEnergySets[k].second = 0;
        }

        // genome 20-25, rnd
        for (int k = 20; k < 25; k++) {
            for (int i = 0; i < controlElectrodeNumber; i++) {
                voltageEnergySets[k].first[controlElectrodeIndices[i]] = enhance::random_double(
                    parameterStorage->parameters.at("controlVoltageMin"),
                    parameterStorage->parameters.at("controlVoltageMax"));
            }
            voltageEnergySets[k].second = 0;
        }

        // mutate
        bool mutatedThisGenome;
        for (int k = 0; k < 25; k++) {
            mutatedThisGenome = false;
            for (int i = 0; i < controlElectrodeNumber; i++) {
                if (enhance::random_double(0, 1) > 0.9) {
                    voltageEnergySets[k].first[controlElectrodeIndices[i]] = enhance::random_triangle(
                        parameterStorage->parameters.at("controlVoltageMin"),
                        voltageEnergySets[k].first[controlElectrodeIndices[i]],
                        parameterStorage->parameters.at("controlVoltageMax"));
                    mutatedThisGenome = true;
                }
            }
            if (mutatedThisGenome) {
                voltageEnergySets[k].second = 0;
            }
        }
        // voltageEnergySets[k].second=0 for all changed voltageEnergySetss (can be
        // used to not run unchanged points (top5) second time, but not implemented
        // yet, bc currents are unknown)

        // ------------ run generation ------
        std::cout << "------------------------------ run geneartion " << generation
                  << " ------------------------------" << std::endl;
        for (int k = 0; k < 25; k++) {
            std::cout << "genome: " << k << " voltages:";
            for (int i = 0; i < controlElectrodeNumber; i++) {
                std::cout << " " << controlElectrodeIndices[i] << ": "
                          << voltageEnergySets[k].first[controlElectrodeIndices[i]];
            }
            std::cout << std::endl;

            auto startTime = std::chrono::steady_clock::now();

            std::pair<std::vector<double>, std::vector<double>> result = jobManager.runControlVoltagesSetup(voltageEnergySets[k].first);
            outputCurrents = result.first;
            outputCurrentUncerts = result.second;

            calcOptimizationEnergy();
            voltageEnergySets[k].second = optEnergy;
            saveResults(k);
            dataFile->addData("generation", &generation);

            std::cout << "optEnergy: " << optEnergy << " fitness: (" << fitness
                      << " +- " << fitnessUncert << ") normedDiff: " << normedDiff
                      << std::endl;
            auto endTime = std::chrono::steady_clock::now();
            std::cout << "time per VoltageSetup = "
                      << std::chrono::duration_cast<std::chrono::milliseconds>(
                             endTime - startTime)
                             .count()
                    / 1000.0
                      << " s" << std::endl;

            if (fitness > bestFitness) {
                bestFitness = fitness;
                bestFitnessUncert = fitnessUncert;
            }
        }

        std::cout << "generation " << generation
                  << " done! sorted results: " << std::endl;
        std::sort(voltageEnergySets.begin(), voltageEnergySets.end(),
            genomeComparator);
        for (int k = 0; k < 25; k++) {
            std::cout << "genome: " << k + 1
                      << " optEnergy: " << voltageEnergySets[k].second
                      << " voltages: ";
            for (int i = 0; i < controlElectrodeNumber; i++) {
                std::cout << " " << controlElectrodeIndices[i] << ": "
                          << voltageEnergySets[k].first[controlElectrodeIndices[i]];
            }
            std::cout << std::endl;
        }

        if ((increaseNumber < parameterStorage->parameters.at("maxStepIncreases")) and (bestFitness + bestFitnessUncert * 2) > 1) {
            increaseNumber++;
            parameterStorage->parameters["calcCurrentSteps"] *= 2;
            std::cout << "############ steps increased!! now: "
                      << parameterStorage->parameters.at("calcCurrentSteps")
                      << " #############" << std::endl;
        }
        if (voltageEnergySets[0].second > parameterStorage->parameters.at("convergenceEnergy") | generation * 25 > parameterStorage->parameters.at("maxIterations")) {
            std::cout
                << "-------------------> optimization stopped <-------------------"
                << std::endl;
            std::cout << "generation: " << generation
                      << " optEnergy: " << voltageEnergySets[0].second << std::endl;
            break;
        }
    }

    DEBUG_FUNC_END
}

/*!
  optimize cotrol voltages using basin hopping
  \param rndStart 0: searchForRandomStart() is called to find best start point
  \n 1: voltages given in input file are used \n 2: used in continue mode,
  voltages have been set before
 */
void Optimizer::optimizeBasinHopping(size_t startMode /*= 0*/)
{
    DEBUG_FUNC_START
    // voltageEnergySets positioning: 0: current. 1: last. 2: current basin best.
    // 3: last basin best
    while (voltageEnergySets.size() < 4) {
        voltageEnergySets.push_back(std::pair<std::vector<double>, double>(
            std::vector<double>(electrodeNumber), 0));
    }

    std::cout << "running optimization - basin hopping" << std::endl;

    double basinAccepted = -1;

    // init voltages
    if (startMode == 0) { // voltages given in input
        for (int i = 0; i < electrodeNumber; i++) {
            voltageEnergySets[0].first[i] = parameterStorage->electrodes[i].voltage;
        }

        std::pair<std::vector<double>, std::vector<double>> result = jobManager.runControlVoltagesSetup(voltageEnergySets[0].first);
        outputCurrents = result.first;
        outputCurrentUncerts = result.second;
        calcOptimizationEnergy();
        voltageEnergySets[0].second = optEnergy;
        saveResults(0);
        dataFile->addData("basinAccepted", &basinAccepted);

        std::cout << "optEnergy: " << voltageEnergySets[0].second << " fitness: ("
                  << fitness << " +- " << fitnessUncert
                  << ") normedDiff: " << normedDiff << std::endl;
        voltageEnergySets[2] = voltageEnergySets[0]; // basin best
        voltageEnergySets[3].second = -INFINITY; // last basin best
    } else if (startMode == 1) { // rnd start
        searchForRandomStart();
    } else if (startMode == 2) { // start set before (continue mode)
    }

    // constructive run
    double lastFitness = fitness;
    double lastFitnessUncert = fitnessUncert;
    double lastNormedDiff = normedDiff;

    voltageEnergySets[1] = voltageEnergySets[0];

    int increaseNumber = 0;
    while (optEnergy < parameterStorage->parameters.at("convergenceEnergy") & iteration < parameterStorage->parameters.at("maxIterations")) {
        iteration++;
        // --------------- basin hopp check ------------------------------
        if (iteration - lastIterationIncrease > parameterStorage->parameters.at(
                "basinWaitSteps")) { // no increase in last basinWaitSteps
            std::cout << "-- stucked in basin! --" << std::endl;
            std::cout << "current basin best: " << voltageEnergySets[2].second
                      << " last basin best: " << voltageEnergySets[3].second
                      << std::endl;

            // check metropolis for last basin
            if ((voltageEnergySets[2].second < voltageEnergySets[3].second) & (enhance::fastExp((voltageEnergySets[2].second - voltageEnergySets[3].second) / parameterStorage->parameters.at("basinTemp")) < enhance::random_double(0, 1))) {
                std::cout << "-- basin not accepted --" << std::endl;
                voltageEnergySets[0] = voltageEnergySets[3]; // update current point (for next basin)
                basinAccepted = 2;
            } else {
                std::cout << "-- basin accepted --" << std::endl;
                voltageEnergySets[0] = voltageEnergySets[2]; // update current point (for next basin)
                voltageEnergySets[3] = voltageEnergySets[2]; // update last basin best
                basinAccepted = 3;
            }

            std::cout << "new random voltages: " << std::endl;
            for (int i = 0; i < controlElectrodeNumber; i++) {
                voltageEnergySets[0]
                    .first[controlElectrodeIndices[i]]
                    = enhance::random_double(
                        std::max(parameterStorage->parameters.at("controlVoltageMin"),
                            voltageEnergySets[0].first[controlElectrodeIndices[i]] - parameterStorage->parameters.at("basinDeltaV")),
                        std::min(parameterStorage->parameters.at("controlVoltageMax"),
                            voltageEnergySets[0].first[controlElectrodeIndices[i]] + parameterStorage->parameters.at("basinDeltaV")));
                std::cout << controlElectrodeIndices[i] << " "
                          << voltageEnergySets[0].first[controlElectrodeIndices[i]]
                          << std::endl;
            }

            // run first point in new basin
            std::pair<std::vector<double>, std::vector<double>> result = jobManager.runControlVoltagesSetup(voltageEnergySets[0].first);
            outputCurrents = result.first;
            outputCurrentUncerts = result.second;
            calcOptimizationEnergy();
            voltageEnergySets[0].second = optEnergy;
            saveResults(0);
            dataFile->addData("basinAccepted", &basinAccepted);

            std::cout << "optEnergy: " << voltageEnergySets[0].second << " fitness: ("
                      << fitness << " +- " << fitnessUncert
                      << ") normedDiff: " << normedDiff << std::endl;

            // seupt for next iteration
            lastFitness = fitness;
            lastFitnessUncert = fitnessUncert;
            lastNormedDiff = normedDiff;
            voltageEnergySets[1] = voltageEnergySets[0];

            voltageEnergySets[2] = voltageEnergySets[0]; // update basin best
            lastIterationIncrease = iteration;

            iteration++;
        }

        // --------------- standard MC ------------------------------
        auto startTime = std::chrono::steady_clock::now();
        // get new random voltages
        std::cout << "new random voltages: " << std::endl;
        for (int i = 0; i < controlElectrodeNumber; i++) {
            voltageEnergySets[0].first[controlElectrodeIndices[i]] = enhance::random_double(
                std::max(parameterStorage->parameters.at("controlVoltageMin"),
                    voltageEnergySets[0].first[controlElectrodeIndices[i]] - parameterStorage->parameters.at("maxDeltaV")),
                std::min(parameterStorage->parameters.at("controlVoltageMax"),
                    voltageEnergySets[0].first[controlElectrodeIndices[i]] + parameterStorage->parameters.at("maxDeltaV")));
            std::cout << controlElectrodeIndices[i] << " "
                      << voltageEnergySets[0].first[controlElectrodeIndices[i]]
                      << std::endl;
        }

        // run
        std::pair<std::vector<double>, std::vector<double>> result = jobManager.runControlVoltagesSetup(voltageEnergySets[0].first);
        outputCurrents = result.first;
        outputCurrentUncerts = result.second;
        calcOptimizationEnergy();
        voltageEnergySets[0].second = optEnergy;
        saveResults(0);

        std::cout << "iteration " << iteration << std::endl;
        std::cout << "now:  optEnergy: " << voltageEnergySets[0].second
                  << " fitness: (" << fitness << " +- " << fitnessUncert
                  << ") normedDiff: " << normedDiff
                  << "\nlast: optEnergy: " << voltageEnergySets[1].second
                  << " fitness: (" << lastFitness << " +- " << lastFitnessUncert
                  << ") normedDiff: " << lastNormedDiff << std::endl;
        if ((voltageEnergySets[0].second < voltageEnergySets[1].second) & (enhance::fastExp((voltageEnergySets[0].second - voltageEnergySets[1].second) / parameterStorage->parameters.at("MCTemp")) < enhance::random_double(0, 1))) {
            std::cout << "-- not accepted --" << std::endl;
            basinAccepted = 0;
            // swap back
            voltageEnergySets[0].first = voltageEnergySets[1].first;
        } else {
            std::cout << "-- accepted --" << std::endl;
            basinAccepted = 1;
            // setup for next iteration
            voltageEnergySets[1] = voltageEnergySets[0];

            lastFitness = fitness;
            lastFitnessUncert = fitnessUncert;
            lastNormedDiff = normedDiff;
        }
        dataFile->addData("basinAccepted", &basinAccepted);

        // check increase
        if (voltageEnergySets[0].second > voltageEnergySets[2].second) {
            voltageEnergySets[2] = voltageEnergySets[0]; // update basin best
            lastIterationIncrease = iteration;
        }

        auto endTime = std::chrono::steady_clock::now();
        std::cout << "time per VoltageSetup = "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(
                         endTime - startTime)
                         .count()
                / 1000.0
                  << " s" << std::endl;

        if ((increaseNumber < parameterStorage->parameters.at("maxStepIncreases")) and (fitness + fitnessUncert * 2) > 1) {
            increaseNumber++;
            parameterStorage->parameters["calcCurrentSteps"] *= 2;
            std::cout << "############ steps increased!! now: "
                      << parameterStorage->parameters["calcCurrentSteps"]
                      << " #############" << std::endl;
        }
    }
    std::cout << "-------------------> optimization stopped <-------------------"
              << std::endl;
    std::cout << "iteration: " << iteration << " optEnergy: " << optEnergy
              << std::endl;
    DEBUG_FUNC_END
}

/*!
  optimize cotrol voltages along the local gradient - not implemented yet
 */
void Optimizer::optimizeGradient()
{
    DEBUG_FUNC_START

    throw std::logic_error("optimizeGradient - not implemented yet");

    // std::vector<double> gradient(controlElectrodeNumber);

    // std::vector<double> voltageEnergySets[1].first =
    // voltageEnergySets[0].first; double voltageEnergySets[1].second = optEnergy;
    // int gradientComponentSign = 1; //if V+deltaV to cals gradient exceeds
    // controlVoltageMax,  V-deltaV ist used instead. info ist stored in this
    // variable double stepWidth = 1; while (true){
    //     //calc Gradient
    //     for (size_t i = 0; i < controlElectrodeNumber; i++){
    //         voltageEnergySets[0].first = voltageEnergySets[1].first;
    //         if (voltageEnergySets[0].first[controlElectrodeIndices[i]] +
    //         parameterStorage->parameters.at("gradDeltaV") <
    //         parameterStorage->parameters.at("controlVoltageMax")){
    //             gradientComponentSign = 1;
    //             voltageEnergySets[0].first[controlElectrodeIndices[i]] +=
    //             parameterStorage->parameters.at("gradDeltaV");
    //         }
    //         else{
    //             gradientComponentSign = -1;
    //             voltageEnergySets[0].first[controlElectrodeIndices[i]] -=
    //             parameterStorage->parameters.at("gradDeltaV");
    //         }

    //         std::pair<std::vector<double>,std::vector<double>> result =
    //         jobManager.runControlVoltagesSetup(voltageEnergySets[0].first);
    //         outputCurrents       = result.first;
    //         outputCurrentUncerts = result.second;

    //         calcOptimizationEnergy();

    //         gradient[i] =
    //         gradientComponentSign*(voltageEnergySets[1].second-optEnergy)/parameterStorage->parameters.at("gradDeltaV");
    //     }

    //     //move Grad Step

    //     if (false){ //converged
    //         break;
    //     }
    // }

    DEBUG_FUNC_END
}

/*!
    generates "rndStartPoints" random control voltage points and sets
   voltageEnergySets[0] to the best result
*/
void Optimizer::searchForRandomStart()
{
    DEBUG_FUNC_START
    // voltageEnergySets positioning: "rndStartPoints" rnd start candidates
    while (voltageEnergySets.size() < parameterStorage->parameters.at("rndStartPoints")) {
        voltageEnergySets.push_back(std::pair<std::vector<double>, double>(
            std::vector<double>(electrodeNumber), 0));
    }

    std::cout << "------ searching for start point ------" << std::endl;

    std::vector<int> controlElectrodeIndices;
    for (int i = 0; i < electrodeNumber; i++) {
        if ((i != parameterStorage->parameters.at("outputElectrode")) & (i != parameterStorage->parameters.at("inputElectrode1")) & (i != parameterStorage->parameters.at("inputElectrode2"))) {
            controlElectrodeIndices.push_back(i);
        }
    }

    auto rndStartPointsComparator =
        [](const std::pair<std::vector<double>, double>& l,
            const std::pair<std::vector<double>, double>& r) {
            return l.second > r.second;
        }; // needed to find best start point

    for (int k = 0; k < parameterStorage->parameters.at("rndStartPoints"); k++) {
        for (int j = 0; j < controlElectrodeNumber; j++) {
            voltageEnergySets[k].first[controlElectrodeIndices[j]] = enhance::random_double(
                parameterStorage->parameters.at("controlVoltageMin"),
                parameterStorage->parameters.at("controlVoltageMax"));
        }
    }

    // run rnd start candidates
    for (int k = 0; k < parameterStorage->parameters.at("rndStartPoints"); k++) {
        std::cout << "rnd start point: " << k << " voltages:";
        for (int i = 0; i < controlElectrodeNumber; i++) {
            std::cout << " " << controlElectrodeIndices[i] << ": "
                      << voltageEnergySets[k].first[controlElectrodeIndices[i]];
        }
        std::cout << std::endl;

        auto startTime = std::chrono::steady_clock::now();

        std::pair<std::vector<double>, std::vector<double>> result = jobManager.runControlVoltagesSetup(voltageEnergySets[k].first);
        outputCurrents = result.first;
        outputCurrentUncerts = result.second;

        calcOptimizationEnergy();
        voltageEnergySets[k].second = optEnergy;
        saveResults(k);
        double a = -1;
        if (optimizationMode == "MC") {
            dataFile->addData("accepted", &a);
        } else if (optimizationMode == "basinHop") {
            dataFile->addData("basinAccepted", &a);
        }

        std::cout << "optEnergy: " << voltageEnergySets[k].second << " fitness: ("
                  << fitness << " +- " << fitnessUncert
                  << ") normedDiff: " << normedDiff << std::endl;
        auto endTime = std::chrono::steady_clock::now();
        std::cout << "time per VoltageSetup = "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(
                         endTime - startTime)
                         .count()
                / 1000.0
                  << " s" << std::endl;
    }

    std::cout << "start search done! sorted results: " << std::endl;
    std::sort(voltageEnergySets.begin(), voltageEnergySets.end(),
        rndStartPointsComparator);
    for (int k = 0; k < parameterStorage->parameters.at("rndStartPoints"); k++) {
        std::cout << "startPoint: " << k + 1
                  << " optEnergy: " << voltageEnergySets[k].second << " voltages: ";
        for (int i = 0; i < controlElectrodeNumber; i++) {
            std::cout << " " << controlElectrodeIndices[i] << ": "
                      << voltageEnergySets[k].first[controlElectrodeIndices[i]];
        }
        std::cout << std::endl;
    }

    fitness = 0;
    fitnessUncert = 0;

    DEBUG_FUNC_END
}

bool Optimizer::desiredLogicFunction(double val1, double val2,
    std::string gate)
{
    DEBUG_FUNC_START

    bool b1 = val1 > parameterStorage->parameters.at("seperationVoltage");
    bool b2 = val2 > parameterStorage->parameters.at("seperationVoltage");

    if (gate == "AND") {
        return (b1 & b2);
    } else if (gate == "NAND") {
        return !(b1 & b2);
    } else if (gate == "OR") {
        return (b1 | b2);
    } else if (gate == "NOR") {
        return !(b1 | b2);
    } else if (gate == "XOR") {
        return (b1 ^ b2);
    } else if (gate == "NXOR") {
        return !(b1 ^ b2);
    } else {
        throw std::runtime_error("logic gate not found");
    }

    DEBUG_FUNC_END
}
