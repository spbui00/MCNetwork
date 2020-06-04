#ifndef JOBHANDLING
#define JOBHANDLING
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
    one job consits of one set of fixed voltages (incl. input voltages). steps to run are split up in tasksPerJob packs. each task can be handeled by single thread.
 */
class Job {
public:
    int ID;
    std::unique_ptr<std::mutex> jobMutex;
    mfem::GridFunction potential; /*!< solution of laplace eq. calculated only once by first thread and saved to pass it to other threads that might join work on this job */
    std::vector<bool> equilOccupation; /*!< equilibrium occupation of system. calculated only once by first thread and saved to pass it to other threads that might join work on this job */

    int equilSteps, totalSteps, stepsPerTask;
    int tasksToGo;
    static const int tasksPerJob = 100;
    int threadNumber = 0; /*!< number of threads that are currently working on this job */
    std::vector<double> voltages;
    double resultCurrent = 0, resultCurrentUncert = 0;

    Job(int ID)
        : ID(ID)
    {
        jobMutex = std::make_unique<std::mutex>();
    };
};

/*!
    handler for parallelization.
 */
class JobManager {
public:
    JobManager(std::shared_ptr<ParameterStorage>);

    std::pair<std::vector<double>, std::vector<double>> const runControlVoltagesSetup(std::vector<double> const& voltages); /*!< returns curr, currUncert */

private:
    int voltageScanPoints, electrodeNumber;
    std::shared_ptr<ParameterStorage> parameterStorage;
    std::vector<System*> systems;
    std::vector<Job> jobs;
    std::mutex jobSearchMutex; /*!< its only allowed to one thread at a time to search for new work! */
    static void handleJobList(std::vector<Job>& jobs, System* const system, std::mutex& searchMutex);
};

#endif // JOBHANDLING_H
