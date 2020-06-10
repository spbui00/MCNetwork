#ifndef SYSTEM_H
#define SYSTEM_H
#include <fstream>
#include <iostream>
#include <memory>
#include <unordered_map>

#include "enhance.hpp"
#include "finiteElemente.h"
#include "parameterstorage.h"

#include "datafile.h"
#include <boost/multi_array.hpp>

#include <shared_mutex>
#include <thread>

#include <chrono>
#include <ctime>

/*!
    main class for physical model
 */
class System {
public:
    int* outputCurrentCounter; /*!< counts hops to/from output electrode */
    double time = 0; /*!< system time. see System::increaseTime() */

    System(std::shared_ptr<ParameterStorage> const&);
    System(System const& oldSys);

    void createRandomNewDevice();
    void loadDevice();

    void initilizeMatrices();
    void getReadyForRun();

    void reset();
    void resetStoredStates();
    void run(int steps);

    void updatePotential(std::vector<double> const& voltages);
    void updatePotential(mfem::GridFunction const& potential);
    void updateOccupationAndPotential(std::vector<bool> const& newOccupation, mfem::GridFunction const& potential);
    std::vector<bool> const& getOccupation();

    mfem::GridFunction getPotential() const;

private:
    int acceptorNumber; /*!< acceptors = dopants*/
    int electrodeNumber; /*!< */
    int hoppingSiteNumber; /*!< acceptors and electrodes are both hoppingSites => hoppingSiteNumber = acceptorNumber + electrodeNumber*/
    double *donorPositionsX, *donorPositionsY, *acceptorPositionsX, *acceptorPositionsY, *electrodePositionsX, *electrodePositionsY; /*!<  2D-array: hopping site distances  1D - array */
    double* distances; /*!<  2D-array: hopping site distances */
    double* energies; /*!<  1D-array: single hopping site energies  */
    int* currentCounter; /*!<  1D-array: counter for all electrodes */
    double* pairEnergies; /*!<  2D-array: coulomb interaction between pair of hopping sites  */
    double* deltaEnergies; /*!<  2D-array: difference between hopping site energies */
    double* rates; /*!<  2D-array: transition rates  */
    double* baseRates; /*!<  2D-array: constant distance dependend part of transition rates (exp(-2r/a)) */
    std::vector<bool> occupation; /*!<  1D-array: occupation of hopping sites */
    std::vector<double> randomEnergies; /*!<  normally distributed site energies caused by various real-world effects */

    std::vector<std::vector<int>> interactionPartners; /*!< indices of coulomb interacting partners. needed bc of cut on coulomb interaction by "maxInteractionDist". sorted by hoppingSite index*/
    std::vector<std::vector<int>> hoppingPartnersAcceptors; /*!< indices of possible hopping partner (acceptors).  needed bc of cut on hoppingPartners by "minHoppingDist" and "maxHoppingDist". sorted by hoppingSite index. example: hoppingPartnersAcceptors[hoppingSiteIndex] <- vector of indices of acceptors where hole form hoppingSiteIndex can hop to. */
    std::vector<std::vector<int>> hoppingPartnersElectrodes; /*!< indices of possible hopping partner (electrodes).  needed bc of cut on hoppingPartners by "minHoppingDist" and "maxHoppingDist". sorted by hoppingSite index. example: hoppingPartnersElectrodes[hoppingSiteIndex] <- vector of indices of electrodes where hole form hoppingSiteIndex can hop to. */

    int lastSwapped1 = 0; /*!< swap 1->2; int = index */
    int lastSwapped2 = 0; /*!< swap 1->2; int = index */

    double ratesSum = 0;
    double constantRatesSumPart = 0; /*!< rates of electrode-electrode hops are constant */
    double locLenA; /*!< physical parameter "a". copied from input file for better performance (otherwise map.at() had to be called a lot of times) */

    std::unique_ptr<FiniteElementeBase> finEle; /*!< finEle device */

    std::shared_ptr<std::vector<double>> partRatesSumList; /*!<list of accumulated rates. only used in storing mode for binary search */
    std::shared_ptr<std::unordered_map<std::vector<bool>, std::shared_ptr<std::vector<double>>>> knownPartRatesSumList; /*!< map of lists of accumulated rates for binary search, to store known states. only used in storing mode. */
    std::shared_ptr<std::unordered_map<std::vector<bool>, double>> knownRatesSum; /*!< map of rate sums for binary search, to store known states. only used in storing mode. */

    std::shared_ptr<ParameterStorage> parameterStorage;
    std::shared_ptr<DataFile> additionalDatafile; /*!< save data in verbose mode*/

    void updateAfterSwap();

    bool readyForRun = false; /*!<set by System::getReadyForRun()  */
    bool storingMode; /*!< if set true performance is optimized by storing known states */
    bool ratesInMemory = false; /*!< save if last step was found in stored states. if true, binary search is done to find swap */

    bool storeKnownStates = true; /*!<set false after memory limit is reached */

    void setOccupation(std::vector<bool> const& newOccupation);

    void findSwap();
    void findSwapBS(); /*!< using binary search */
    void updateRatesStoringMode(); /*!< single processor, storing mode */
    void updateRates(); /*!< core calculation of rates */
    void increaseTime();

    void resetPotential();
    void setNewPotential();
};

#endif // SYSTEM_H
