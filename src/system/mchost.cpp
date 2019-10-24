#include "mchost.h"
#ifndef NDEBUG
   #include "debug.h"
#endif


MCHost::MCHost(std::shared_ptr<ParameterStorage> _parameterStorage)
{
    #ifndef NDEBUG
        DEBUG_FUNC_START
    #endif

    parameterStorage=_parameterStorage;

    #ifndef NDEBUG
        DEBUG_FUNC_END
    #endif
}

void MCHost::setup()
{
    #ifndef NDEBUG
        DEBUG_FUNC_START
    #endif
    
  
    #ifndef NDEBUG
        DEBUG_FUNC_END
    #endif
}

