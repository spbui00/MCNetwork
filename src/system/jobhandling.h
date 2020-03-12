#ifndef JOBHANDLING
#define JOBHANDLING
#include <iostream>
#include <memory>
#include <algorithm>
#include <fstream>
#include "parameterstorage.h"
#include "system.h"
#include <float.h>
#include <chrono>
#include <vector>
#include <map>

//(un)comment to (en/dis)able time tracker
#define TIMETRACKER


class Job
{
    public:
        int ID; 
        std::unique_ptr<std::mutex> jobMutex;
        mfem::GridFunction potential;
        std::vector<bool> equilOccupation;

        int equilSteps, totalSteps, stepsPerTask;
        int tasksToGo; 
        static const int tasksPerJob = 20; //number of 
        int threadNumber = 0;
        std::vector<double> voltages;
        double resultCurrent = 0, resultCurrentUncert = 0;

        Job(int ID): ID(ID) {jobMutex = std::make_unique<std::mutex>();};

        #ifdef TIMETRACKER
            std::ofstream timeFile;
        #endif
};


class JobManager
{
public:
    JobManager(std::shared_ptr<ParameterStorage>);

    std::pair<std::vector<double>,std::vector<double>> const runControlVoltagesSetup(std::vector<double> const & voltages); //returns curr, currUncert

private:
    int voltageScanPoints, electrodeNumber;
    std::shared_ptr<ParameterStorage> parameterStorage;
    std::vector<System * > systems;
    std::vector<Job> jobs;
    std::mutex jobSearchMutex;
    static void handleJobList(std::vector<Job> & jobs, System * const system, std::mutex & searchMutex);

};





#endif // JOBHANDLING_H
