#ifndef SYSTEM_H
#define SYSTEM_H
#include <iostream>
#include <memory>
#include <map>
#include <fstream>
// #include "datafile.h"
#include "parameterstorage.h"
#include <chrono>
#include <ctime>

class System
{
private:
    int steps=0;

    std::shared_ptr<ParameterStorage> parameterStorage;
    // std::unique_ptr<DataFile> dataFile;
    
    
public:
    void setup();
    System(std::shared_ptr<ParameterStorage>);

};

#endif // SYSTEM_H
