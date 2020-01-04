#include "../lib/enhance.hpp"

#include "system/mchost.h"

#include "debug.h"

#include <random>
#include <iostream>
#include <boost/program_options.hpp>


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
      ("mnd", "make new device")
      // ("dir", po::value<std::string>(),"define working dir. has to contain 'in.txt'")s
      ("help", "produce help message");
   
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
      workingDirecotry ="../data/";
   }
   std::cout<<"running in dir: "<<workingDirecotry<<std::endl;
   
   std::string inputFileName=workingDirecotry + "in.txt";
   std::shared_ptr<ParameterStorage> parameterStorage(new ParameterStorage(inputFileName));//all input parameters are stored in the shared pointer "inputfile". all classes get the pointer 
   parameterStorage->workingDirecotry=workingDirecotry;

   MCHost mchost(parameterStorage);
   mchost.setup(vm.count("mnd"));
   mchost.optimizeMC();

   DEBUG_FUNC_END
   return 0;
}
