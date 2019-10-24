#include "../lib/enhance.hpp"
#include "system/mchost.h"

#ifndef NDEBUG
   #include "debug.h"
#endif

#include <random>
#include <iostream>
#include <boost/program_options.hpp>




int main(int argc, char *argv[]){
   #ifndef NDEBUG
      DEBUG_FUNC_START
   #endif

   std::string inputFileName="in.txt";
   std::shared_ptr<ParameterStorage> parameterStorage;
   parameterStorage.reset(new ParameterStorage(inputFileName));  //all input parameters are stored in the shared pointer "inputfile". all classes get the pointer 

   std::cout<<parameterStorage->parameters.at("T")<<std::endl;


   MCHost mchost(parameterStorage);
   
   #ifndef NDEBUG
      DEBUG_FUNC_END
   #endif

   return 0;
}
