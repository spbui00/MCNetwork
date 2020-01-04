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

class MCHost
{
private:
    int steps=0;
    int hoppingSiteNumber=0;
    int voltageScanPointsNumber;
    int electrodeNumber;
    double fitness;
    double ratesSum=0;
    double locLenA;
    double** rates;
    double* outputCurrentBuffer;


    std::shared_ptr<ParameterStorage> parameterStorage;
    std::shared_ptr<DataFile>         dataFile;
    std::unique_ptr<System>           system;

    void makeSwap();
    void calcRates();
    void calcFitness();
    void saveResults();
    bool desiredLogicFunction(double val1, double val2, std::string gate);

    
public:
    MCHost(std::shared_ptr<ParameterStorage>);
    void setup(bool makeNewDevice=true);

    void singleRun();
    void runVoltageSetup();
    void optimizeMC();
    void run();
    

};

#endif // MCHOST_H
