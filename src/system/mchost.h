#ifndef MCHOST_H
#define MCHOST_H
#include <iostream>
#include <memory>
#include <fstream>
#include "parameterstorage.h"
#include "system.h"
#include "datafile.h"
#include <float.h>

class MCHost
{
private:
    std::string const mode = "AND";
    int steps=0;
    int hoppingSiteNumber=0;
    int voltageScanPointsNumber;
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

    
public:
    MCHost(std::shared_ptr<ParameterStorage>);
    void setup(bool makeNewDevice=true);

    void singleRun();
    void runControlSetup();
    void run();

};

#endif // MCHOST_H
