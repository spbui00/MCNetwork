#ifndef PARAMETERSTORAGE_H
#define PARAMETERSTORAGE_H
#include<fstream>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include <cmath>
#include "debug.h"


struct ElectrodeParameters{
    double pos;
    int edge;
    double voltage;
};


class ParameterStorage
{
    
public:
    std::string gate;
    std::map<std::string,double> parameters; //general parameter map
    std::vector<ElectrodeParameters> electrodes;
    std::string workingDirecotry ="./";

    ParameterStorage(std::string);
};



#endif // PARAMETERSTORAGE_H
