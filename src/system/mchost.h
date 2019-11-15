#ifndef MCHOST_H
#define MCHOST_H
#include <iostream>
#include <memory>
#include <fstream>
#include "parameterstorage.h"
#include "system.h"

class MCHost
{
private:
    int steps=0;
    int acceptorNumber=0;
    double ratesSum=0;
    double locLenA=0;
    double** rates;

    std::shared_ptr<ParameterStorage> parameterStorage;
    std::unique_ptr<System> system;

    void makeSwap();
    void calcRates();

    
public:
    MCHost(std::shared_ptr<ParameterStorage>);
    void setup();
    void run(int N);
};

#endif // MCHOST_H
