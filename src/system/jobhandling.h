#ifndef JOBHANDLING
#define JOBHANDLING
#include <iostream>
#include <memory>
#include <algorithm>
#include <fstream>
#include "parameterstorage.h"
#include "system.h"
#include "datafile.h"
#include <float.h>
#include <chrono>
#include <vector>
#include <map>


class Job
{
    public:
        // Job(){};
        int equilSteps, totalSteps, stepsPerTask;
        int tasksToGo; 
        static const int tasksPerJob = 20; //number of 
        int threadNumber = 0;
        std::vector<double> voltages;
        double resultCurrent = 0, resultCurrentUncert = 0;
};


class JobManager
{
private:
    int voltageScanPoints, electrodeNumber;
    std::shared_ptr<ParameterStorage> parameterStorage;
    std::vector<System * > systems;
    std::vector<Job> jobs;
    std::mutex jobSearchMutex;
    static void handleJobList(std::vector<Job> & jobs, System * system, std::mutex & searchMutex);


    void singleRun();
public:
    JobManager(std::shared_ptr<ParameterStorage>);

    std::pair<std::vector<double>,std::vector<double>> runControlVoltagesSetup(std::vector<double> voltages); //returns curr, currUncert
};





#endif // JOBHANDLING_H
