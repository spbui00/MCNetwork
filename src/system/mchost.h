#ifndef MCHOST_H
#define MCHOST_H
#include <iostream>
#include <memory>
#include <fstream>
#include "parameterstorage.h"

class MCHost
{
private:
    int steps=0;

    std::shared_ptr<ParameterStorage> parameterStorage;
    
    
public:
    MCHost(std::shared_ptr<ParameterStorage>);
    void setup();
};

#endif // MCHOST_H
