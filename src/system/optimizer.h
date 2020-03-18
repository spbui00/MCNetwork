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
    void optimizeGenetic(std::vector<std::pair<std::vector<double>,double>> const & startGenome = {});
    void optimizeBasinHopping(bool rndStart = false);
    void optimizeGradient(int basinNumber);
    void run();
    
    
private:
    int electrodeNumber, voltageScanPoints;
    int controlElectrodeNumber;
    double fitness = 0,fitnessUncert = 0,optEnergy = 0,normedDiff = 0;

    std::string mode;

    std::vector<std::vector<double>> voltageSets;
    std::vector<std::vector<double>> outputCurrents;
    std::vector<std::vector<double>> outputCurrentUncerts;
    std::vector<int> controlElectrodeIndices;
    

    std::shared_ptr<ParameterStorage> parameterStorage;
    std::shared_ptr<DataFile>         dataFile;
    
    std::shared_ptr<std::vector<System * >> systems;
    JobManager jobManager;

    void calcOptimizationEnergy();
    void saveResults();
    bool desiredLogicFunction(double val1, double val2, std::string gate);

    void searchForRandomStart();
};


#endif // OPTIMIZER_H
