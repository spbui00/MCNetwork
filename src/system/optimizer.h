#ifndef OPTIMIZER
#define OPTIMIZER
#include "datafile.h"
#include "jobhandling.h"
#include "parameterstorage.h"
#include "system.h"
#include <algorithm>
#include <chrono>
#include <float.h>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

/*!
    container for optimization routines.
 */
class Optimizer {
public:
    Optimizer(std::shared_ptr<ParameterStorage>);

    void run(std::string optimizationMode, int startMode);

private:
    int electrodeNumber; /*!< number of electrodes */
    int voltageScanPoints; /*!< number of inout voltages apllied to one input electrode. i.e. 2 for 2x2 search */
    size_t controlElectrodeNumber; /*!< number of control electrodes = electrodeNumber - 3 */
    double fitness = 0; /*!< raw fitness */
    double fitnessUncert = 0; /*!< uncertainty of raw fitness*/
    double optEnergy = 0; /*!< !! optimization energy is called "(corrected) fitness/ \mathcal{F}" in thesis!!. optimization energy of last simulated voltage set. set by calcOptimizationEnergy(). only local variable, use voltageEnergySets[i].second instead. */
    double normedDiff = 0; /*!< parameter to quantify difference between high und low output current */
    size_t iteration = 0;
    size_t lastIterationIncrease = 0; /*!< only used in basinHop mode. last time optEnergy increased */

    std::string optimizationMode; /*!< MC, genetic, basinHop, singleRun */

    std::vector<std::pair<std::vector<double>, double>> voltageEnergySets; /*!< all methods share one central storage for all sets of control voltages and corresponding optEnergies they need to know at a time. size and indexing differs from method to method, see methods doc for more details.\n voltageEnergySets consists of pair of vector of control voltages and optEnergy. \n first = voltages \n second = optEnergy. */
    std::vector<double> outputCurrents; /*!< outputCurrents for one set of control electrodes */
    std::vector<double> outputCurrentUncerts; /*!< outputCurrent uncertainties for one set of control electrodes */

    std::vector<size_t> controlElectrodeIndices; /*!< [i for i in range(electrodeNumber) if (i != outputElectrode and i != inputElectrode)] */

    std::shared_ptr<ParameterStorage> parameterStorage;
    std::shared_ptr<DataFile> dataFile;

    std::shared_ptr<std::vector<System*>> systems; /*!< multiple systems are created for parallelization and stored in shared pointer to be accessiable by jobManeger and optimizer */
    JobManager jobManager;

    void calcOptimizationEnergy();
    void saveResults(size_t index = 0);
    bool desiredLogicFunction(double val1, double val2, std::string gate);

    void searchForRandomStart();

    void singleRun(size_t startMode = 0);
    void optimizeMC(size_t startMode = 0);
    void optimizeGenetic(size_t startMode = 0);
    void optimizeBasinHopping(size_t startMode = 0);
    void optimizeGradient();

    void continueSimulation();
};

#endif // OPTIMIZER_H
