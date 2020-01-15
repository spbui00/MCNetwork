#ifndef MCHOST_H
#define MCHOST_H
#include <iostream>
#include <memory>
#include <algorithm>
#include <fstream>
#include "parameterstorage.h"
#include "system.h"
#include "datafile.h"
#include <float.h>
#include <chrono>
#include <vector>
#include <map>

class MCHost
{
private:
    int hoppingSiteNumber=0;
    int voltageScanPointsNumber;
    int electrodeNumber;
    double fitness,fitnessUncert,optEnergy,normedDiff;

    double outputCurrent,outputCurrentSqrt,outputCurrentStd;
    double * outputCurrentBuffer, * outputCurrentUncertBuffer;

    // std::map<std::string,std::shared_ptr<std::vector<double>>> knownRates;
    // std::map<std::string,double>  knownRatesSum;



    std::shared_ptr<ParameterStorage> parameterStorage;
    std::shared_ptr<DataFile>         dataFile;
    std::unique_ptr<System>           system;

    void calcOptimizationEnergy();
    void saveResults();
    bool desiredLogicFunction(double val1, double val2, std::string gate);


public:
    MCHost(std::shared_ptr<ParameterStorage>);
    void setup(bool makeNewDevice=true);

    void singleRun();
    void runVoltageSetup();
    void optimizeMC(bool rndStart = false);
    void run();
    

};

#endif // MCHOST_H
