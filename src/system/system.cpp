#include "system.h"
#ifndef NDEBUG
   #include "debug.h"
#endif


System::System(std::shared_ptr<ParameterStorage> _parameterStorage)
{
    #ifndef NDEBUG
        DEBUG_FUNC_START
    #endif

    parameterStorage=_parameterStorage;

    #ifndef NDEBUG
        DEBUG_FUNC_END
    #endif
}

void System::setup()
{
    #ifndef NDEBUG
        DEBUG_FUNC_START
    #endif

  
    #ifndef NDEBUG
        DEBUG_FUNC_END
    #endif
}

