#include "../lib/enhance.hpp"

#include "system/mchost.h"

#include "debug.h"

#include <random>
#include <iostream>
#include <boost/program_options.hpp>




int main(int argc, char *argv[]){
   DEBUG_FUNC_START

   enhance::seed = std::random_device{}();
   #ifndef NDEBUG
      enhance::seed = 123456749;
   #endif
   enhance::rand_engine.seed(enhance::seed);

   std::string workingDirecotry ="../data/";

   std::string inputFileName=workingDirecotry + "in.txt";
   std::shared_ptr<ParameterStorage> parameterStorage(new ParameterStorage(inputFileName));//all input parameters are stored in the shared pointer "inputfile". all classes get the pointer 
   parameterStorage->workingDirecotry=workingDirecotry;

   MCHost mchost(parameterStorage);
   mchost.setup(false);
   mchost.run();


   


   
   DEBUG_FUNC_END
   return 0;
}
