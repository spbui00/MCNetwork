#ifndef SYSTEM_H
#define SYSTEM_H
#include <iostream>
#include <memory>
#include <unordered_map>
#include <fstream>
// #include "datafile.h"
#include "parameterstorage.h"
#include "../lib/enhance.hpp"
#include "../lib/finiteElemente/finiteElemente.h"

#include <boost/multi_array.hpp>



#include <chrono>
#include <ctime>

class System
{
private:
    int acceptorNumber, hoppingSiteNumber, electrodeNumber;
    double * donorPositionsX, * donorPositionsY, * acceptorPositionsX, * acceptorPositionsY, * electrodePositionsX, * electrodePositionsY; // 1D
    double * energies; // 1D
    bool * occupation; // 1D
    double * pairEnergies, * distances, * deltaEnergies; //2D
    unsigned long long hasedCurrentState=0; //state is stored hased (better performance compared to string)
    
    int lastSwapped1,lastSwapped2; // swap 1->2 int = index

    double ratesSum=0;
    double locLenA;
    std::shared_ptr<std::vector<double>> rates;

    std::shared_ptr<ParameterStorage> parameterStorage;
    
    FiniteElemente * finEle; //finEle device

    void updateAfterSwap();
    
public:

    std::unordered_map<unsigned long long,std::shared_ptr<std::vector<double>>> knownRates;
    std::unordered_map<unsigned long long,double>  knownRatesSum;

    double * currentCounter; // 1D

    double time=0;

    System(std::shared_ptr<ParameterStorage>);



    void createRandomNewDevice();
    void loadDevice();


    void initilizeMatrices();
    void getReadyForRun();

    void makeSwap();
    void updateRates(bool storeKnowStates = false);
    void updatePotential();
    void increaseTime();
};




#endif // SYSTEM_H
