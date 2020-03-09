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
    double pos; //used as angle in circle mode, or as length in rect mode
    int edge;   //only used in circle mode
    double voltage;
};


class ParameterStorage
{
    
public:
    std::string gate;
    std::string geometry;
    std::map<std::string,double> parameters; //general parameter map
    std::vector<ElectrodeParameters> electrodes;
    std::string workingDirecotry ="./";
    
    bool makeNewDevice = false;
    std::vector<double> inputVoltages; //input voltage points to scan


    ParameterStorage(std::string);
};



#endif // PARAMETERSTORAGE_H
