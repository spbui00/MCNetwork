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

    void optimizeMC(bool rndStart = false);
    void optimizeGenetic();
    void run();
    
    
private:
    int electrodeNumber, voltageScanPoints;
    double fitness = 0,fitnessUncert = 0,optEnergy = 0,normedDiff = 0;

    std::vector<std::vector<double>> voltageSets;
    std::vector<std::vector<double>> outputCurrents;
    std::vector<std::vector<double>> outputCurrentUncerts;
    

    std::shared_ptr<ParameterStorage> parameterStorage;
    std::shared_ptr<DataFile>         dataFile;
    
    std::shared_ptr<std::vector<System * >> systems;
    JobManager jobManager;

    void calcOptimizationEnergy();
    void saveResults();
    bool desiredLogicFunction(double val1, double val2, std::string gate);
};


#endif // OPTIMIZER_H
