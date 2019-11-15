#include "dopant.h"
#include "debug.h"


Dopant::Dopant(double posX_, double posY_)
{
    DEBUG_FUNC_START
    
    posX=posX_;
    posY=posY_;


    DEBUG_FUNC_END
}

Dopant::Dopant()
{
    DEBUG_FUNC_START
    
    posX=1;
    posY=1;


    DEBUG_FUNC_END
}

