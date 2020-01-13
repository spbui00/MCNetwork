#ifndef HOPPINGSITE_H
#define HOPPINGSITE_H
#include <iostream>
#include <memory>
#include <fstream>
#include "parameterstorage.h"



class HoppingSite{
//electrodes and dopants are both hopping sites, but only dopants can change their occupation. occupation change of electrodes has no effect
protected:
    bool occupied=false;
public:
    HoppingSite(double posX, double posY): posX(posX), posY(posY){}
    double const posX,posY;
    int currentCounter=0;
    // int absCurrentCounter=0;
    double energy=0,energy2=0,constEnergyPart=0;

    virtual void setOccupation(bool occ){occupied = occ;}
    bool const & getOccupation() const {return occupied;}
};

class Dopant : public HoppingSite{
    using HoppingSite::HoppingSite;
public:
    void setOccupation(bool occ){occupied = occ;}
};

class Electrode : public HoppingSite{
public:
    Electrode(double posX, double posY): HoppingSite(posX, posY) {occupied=true;}
    void setOccupation(bool occ){}
};

#endif // HOPPINGSITE_H
