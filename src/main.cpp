#include "enhance.hpp"

#include "system/optimizer.h"

#include "debug.h"

#include <random>
#include <iostream>
#include <boost/program_options.hpp>

/*! \mainpage MCNetwork
 *
 * \section install_sec Installation
 *
 * \subsection step1 Step 1: Dependencies
 * mfem (serial version): https://mfem.org/building/ \n 
 * boost \n 
 * 
 * \subsection step2 Step 2: Install
 * 
 * git clone https://github.com/MarlonBecker/MCNetwork \n 
 * mkdir build \n 
 * cd build \n 
 * cmake .. \n 
 * make \n 
 * 
 * \subsection step3 Step 3: Run Example
 * 
 * cd data \n 
 * ../build/MCnetwork --mnd --optMC --rSV \n 
 */


/*
TODO:
   TECHNICAL:
      - set "-D SWAPTRACKER" in cmake file
      - improve argument parsing
      - add log file
      - support INT data format in datafile
      - split up optimizer in child classes
      - make one central hdf5 file to combine SWAPTRACKER and TiMETRACKER
   PHYSICAL:
   ANALYSIS:
OPT:
   - count only output electrode current
   - improve hash algorithm for storing mode
   - in optimzer: store only control voltages
   - in genetic algorithm: do not run unchanged genomes again
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
      ("continue"   , "continues last optimization. optimization mode will be detected, no further options needed")
      ("mnd"        , "make new device")
      ("optMC"      , "optimize control voltages using Monte Carlo search")
      ("optGen"     , "optimize control voltages using genetic algorithm")
      ("optBasinHop", "optimize control voltages using basin hopping")
      ("run"        , "just run control voltages defined in in.txt")
      ("rSV"        , "random start voltages (only in combination with opt)")
      ("dir"        , po::value<std::string>(),"define working dir. has to contain 'in.txt'")
      ("verbose"    , "additionally outputs occupation, swapps and time to 'additionalData.hdf5'")
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
   parameterStorage->verbose       = vm.count("verbose");
   int startMode                   = vm.count("rSV");


   if (parameterStorage->verbose){
      if (parameterStorage->parameters.at("threads") > 1) throw std::logic_error("verbose mode is only supported in single processor mode");
      parameterStorage->parameters["additionalFileNumber"] = 0;
   }



   std::string optimizationMode;
   if (vm.count("optMC")){
      optimizationMode = "MC";
   }
   else if (vm.count("optGen")){
      optimizationMode = "genetic";
   }
   else if (vm.count("optBasinHop")){
      optimizationMode = "basinHop";
   }
   else if (vm.count("run")){
      optimizationMode = "singleRun";
   }
   if (vm.count("continue")){
      optimizationMode = "continue";
   }

   Optimizer optimizer(parameterStorage);
   optimizer.run(optimizationMode, startMode);

      

   DEBUG_FUNC_END
   return 0;
}
