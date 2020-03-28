#ifndef OPTIMIZER
#define OPTIMIZER
#include <iostream>
#include <memory>
#include <algorithm>
#include <fstream>
#include "parameterstorage.h"
#include "system.h"
#include "jobhandling.h"
#include "datafile.h"
#include <float.h>
#include <chrono>
#include <vector>
#include <map>

class Optimizer
{
public:
    Optimizer(std::shared_ptr<ParameterStorage>);

    void run(std::string optimizationMode, int startMode);

private:
    int electrodeNumber, voltageScanPoints;
    size_t controlElectrodeNumber;
    double fitness = 0,fitnessUncert = 0,optEnergy = 0;  /*!< optimization energy of last simulated voltage set. set by calcOptimizationEnergy(). only local variable, use voltageEnergySets[i].second instead. */
    double normedDiff = 0; /*!< parameter to quantify difference between high und low output current */
    size_t iteration=0;
    size_t lastIterationIncrease = 0; /*!< only used in basinHop mode. last time optEnergy increased */

    std::string optimizationMode;  /*!< MC, genetic, basinHop, singleRun */

    std::vector<std::pair<std::vector<double>,double>> voltageEnergySets; /*!< all methods share one central storage for all sets of control voltages and corresponding optEnergies they need to know at a time. size and indexing differs from method to method, see methods doc for more details.\n voltageEnergySets consists of pair of vector of control voltages and optEnergy. \n first = voltages \n second = optEnergy. */
    std::vector<double> outputCurrents;
    std::vector<double> outputCurrentUncerts;
    
    std::vector<size_t> controlElectrodeIndices;
    

    std::shared_ptr<ParameterStorage> parameterStorage;
    std::shared_ptr<DataFile>         dataFile;
    
    std::shared_ptr<std::vector<System * >> systems;
    JobManager jobManager;

    void calcOptimizationEnergy();
    void saveResults(size_t index = 0);
    bool desiredLogicFunction(double val1, double val2, std::string gate);

    void searchForRandomStart();

    void optimizeMC          (size_t startMode = 0);
    void optimizeGenetic     (size_t startMode = 0);
    void optimizeBasinHopping(size_t startMode = 0);
    void optimizeGradient();
    void singleRun();

    void continueSimulation();
};


#endif // OPTIMIZER_H
