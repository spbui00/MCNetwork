#ifndef DOPANT_H
#define DOPANT_H
#include <iostream>
#include <memory>
#include <fstream>
#include "parameterstorage.h"

class Dopant
{
private:
    int steps=0;

public:
    Dopant(double posX_, double posY_);
    Dopant();
    double posX,posY;
    double energy=0,constEnergyPart=0;
    bool occupied=false;

};

#endif // DOPANT_H
