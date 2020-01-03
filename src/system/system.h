#ifndef SYSTEM_H
#define SYSTEM_H
#include <iostream>
#include <memory>
#include <map>
#include <fstream>
// #include "datafile.h"
#include "parameterstorage.h"
#include "hoppingSite.h"
#include "../lib/enhance.hpp"
#include "../lib/finiteElemente/finiteElemente.h"

#include <boost/multi_array.hpp>



#include <chrono>
#include <ctime>

class System
{
private:
    int steps, acceptorNumber, hoppingSiteNumber;
    double** pairEnergies;
    double** donorPositions;
    std::shared_ptr<ParameterStorage> parameterStorage;
    
    FiniteElemente * finEle; //finEle device
    // std::unique_ptr<DataFile> dataFile;
    
public:
    System(std::shared_ptr<ParameterStorage>);
    double** distances;
    double** deltaEnergies;


    std::vector< HoppingSite * > hoppingSites {};


    void createRandomNewDevice();
    void loadDevice();

    void initilizeMatrices();
    void getReadyForRun();

    void calcEnergies();
    void setElectrodeVoltage(int electrodeIndex, double voltage);
    void updatePotential();
};




#endif // SYSTEM_H
