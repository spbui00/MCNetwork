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
    double pos; /*!< used as angle (in degrees) in circle mode, or as side position (in fractions of side length) in rect mode */
    int edge;   /*!< only used in rect mode: \n 
                0 --> x = 0,          y = pos * lenY \n 
                1 --> x = lenX,       y = pos * lenY \n 
                2 --> x = pos * lenX, y = 0 \n 
                2 --> x = pos * lenX, y = lenY */
    double voltage;
};


class ParameterStorage
{
    
public:
    std::string gate;
    std::string geometry;
    std::map<std::string,double> parameters; /*!< general parameter map */
    std::vector<ElectrodeParameters> electrodes;
    std::string workingDirecotry ="./";
    
    bool makeNewDevice = false;
    std::vector<double> inputVoltages; /*!< input voltage points to scan */


    ParameterStorage(std::string);
};



#endif // PARAMETERSTORAGE_H
