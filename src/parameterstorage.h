#ifndef PARAMETERSTORAGE_H
#define PARAMETERSTORAGE_H
#include<fstream>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include <cmath>
#include "debug.h"


class ParameterStorage
{
    
public:
    double T, kBT;
    std::map<std::string,double> parameters; //general parameter map

    ParameterStorage(std::string);
};

#endif // PARAMETERSTORAGE_H
