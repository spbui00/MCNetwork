#ifndef SYSTEM_H
#define SYSTEM_H
#include <iostream>
#include <memory>
#include <map>
#include <fstream>
// #include "datafile.h"
#include "parameterstorage.h"
#include "dopant.h"
#include "../lib/enhance.hpp"
#include <boost/multi_array.hpp>



#include <chrono>
#include <ctime>

class System
{
private:
    int steps=0;
    int acceptorNumber=0;
    double** pairEnergies;
    // bool* notOccVec;

    std::shared_ptr<ParameterStorage> parameterStorage;
    // std::unique_ptr<DataFile> dataFile;
    
    
public:
    System(std::shared_ptr<ParameterStorage>);
    double** distances;
    double** deltaEnergies;

    std::vector<Dopant> dopants {};

    void setup();
    void calcEnergies();


};

#endif // SYSTEM_H
