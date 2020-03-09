#include "../lib/enhance.hpp"

#include "system/optimizer.h"

#include "debug.h"

#include <random>
#include <iostream>
#include <boost/program_options.hpp>


/*
TODO:
   TECHNICAL:
      - set "-D SWAPTRACKER" in cmake file
      - improve argument parsing
      - make crashed simulaten resumeable
      - add log file
      - support INT data format in datafile
   PHYSICAL:
      - implement replica exchange
      - implement basin hopping
      - implement circular area
   ANALYSIS:
      - calc correlation between voltages and visible currents (using swap track)
      - calc single trajectories
OPT:
   - store only output electrode current
   - improve hash algorithm
   - do not recalc all rates (if coulomb cut)
*/

namespace po = boost::program_options;


int main(int argc, char *argv[]){
   DEBUG_FUNC_START

   enhance::seed = std::random_device{}();
   #ifndef NDEBUG
      enhance::seed = 123456749;
   #endif
   enhance::rand_engine.seed(enhance::seed);




   po::options_description desc("Allowed options");
   desc.add_options()
      ("mnd"        , "make new device")
      ("optMC"      , "optimize control voltages using Monte Carlo Search")
      ("optGen"     , "optimize control voltages using GeneticAlgorithm")
      ("optBasinHop", "optimize control voltages using Basin Hopping")
      ("run"        , "just run control voltages defined in in.txt")
      ("rSV"        , "random start voltages (only in combination with opt)")
      ("dir"        , po::value<std::string>(),"define working dir. has to contain 'in.txt'")
      ("help"       , "produce help message");

   po::variables_map vm;
   po::store(po::parse_command_line(argc, argv, desc), vm);


   if (vm.count("help"))
   {
      std::cout << desc << "\n";
      return 1;
   }

   std::string workingDirecotry; 
   if (vm.count("dir")){
      workingDirecotry= vm["dir"].as<std::string>();
   }
   else{
      workingDirecotry ="./";
   }
   std::cout<<"running in dir: "<<workingDirecotry<<std::endl;

   std::string inputFileName=workingDirecotry + "in.txt";
   std::shared_ptr<ParameterStorage> parameterStorage(new ParameterStorage(inputFileName));//all input parameters are stored in the shared pointer "inputfile". all classes get the pointer 
   parameterStorage->workingDirecotry=workingDirecotry;
   parameterStorage->makeNewDevice = vm.count("mnd");


   Optimizer optimizer(parameterStorage);


   if (vm.count("optMC")){
      optimizer.optimizeMC(vm.count("rSV"));
   }
   else if (vm.count("optGen")){
      optimizer.optimizeGenetic();
   }
   else if (vm.count("optBasinHop")){
      optimizer.optimizeBasinHopping(vm.count("rSV"));
   }
   else if (vm.count("run")){
      optimizer.run();
   }
      

   DEBUG_FUNC_END
   return 0;
}
