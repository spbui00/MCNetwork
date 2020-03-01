#ifndef SYSTEM_H
#define SYSTEM_H
#include <iostream>
#include <memory>
#include <unordered_map>
#include <fstream>

#include "parameterstorage.h"
#include "../lib/enhance.hpp"
#include "../lib/finiteElemente/finiteElemente.h"

#include <boost/multi_array.hpp>

#include <thread> 
#include <shared_mutex>


#include <chrono>
#include <ctime>


//(un)comment to (en/dis)able swap tracker
// #define SWAPTRACKER

class System
{
public:
    int * outputCurrentCounter;
    double time=0;


    System(std::shared_ptr<ParameterStorage> const &);
    System(System const & oldSys);
    
    void createRandomNewDevice();
    void loadDevice();

    void initilizeMatrices();
    void getReadyForRun();

    void reset();
    void resetStoredStates();
    void run(int steps);
    
    void updatePotential(std::vector<double> const & voltages);
    void updatePotential(mfem::GridFunction  const & potential);
    void updateOccupationAndPotential(std::vector<bool> const & newOccupation, mfem::GridFunction  const & potential);
    std::vector<bool> const & getOccupation();


    mfem::GridFunction getPotential() const;


private:
    int acceptorNumber, hoppingSiteNumber, electrodeNumber;
    double * donorPositionsX, * donorPositionsY, * acceptorPositionsX, * acceptorPositionsY, * electrodePositionsX, * electrodePositionsY; // 1D
    double * energies; // 1D
    int * currentCounter; // 1D
    std::vector<bool> occupation; // 1D
    double * pairEnergies, * rates, * distances, * deltaEnergies; //2D
    
    int lastSwapped1=0,lastSwapped2=0; // swap 1->2 int = index

    double ratesSum=0;
    double constantRatesSumPart=0;
    double locLenA;

    std::unique_ptr<FiniteElemente> finEle; //finEle device

    shared_ptr<std::vector<double>> partRatesSumList; //list of accumulated rates for binary search
    std::shared_ptr< std::unordered_map<std::vector<bool>,std::shared_ptr<std::vector<double>>>> konwnPartRatesSumList; //map of lists of accumulated rates for binary search, to store known states
    std::shared_ptr< std::unordered_map<std::vector<bool>,double>>  knownRatesSum;

    std::shared_ptr<ParameterStorage> parameterStorage;
    

    void updateAfterSwap();

    bool readyForRun=false;
    bool storingMode; // if set true performance is optimized by storing known states
    bool ratesInMemory =false; //save if last step was found in stored states. if true, binary search is done to find swap


    bool storeKnownStates = true;


    #ifdef SWAPTRACKER
        ofstream swapTrackFile;
        int fileNumber=1;
    #endif

    void setOccupation(std::vector<bool> const & newOccupation);


    void findSwap();
    void findSwapBS(); //using binary search
    void updateRatesStoringMode(); //single processor, storing mode
    void updateRates();// core calculation of rates
    void increaseTime();

    void resetPotential();
    void setNewPotential();
};




#endif // SYSTEM_H
