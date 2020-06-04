#ifndef NDEBUG
#include "debug.h"

int myDebugger::funcTraceLevel = -1;
std::set<std::string> myDebugger::debugKeys = { "funcTrace" };
// {"funcTrace"}

void myDebugger::debugPrint(std::string message, std::string key)
{
    if (debugKeys.find(key) != debugKeys.end()) {
        if (key == "bla") {
        } else if (key == "blubb") {
        } else {
            std::cout << message << std::endl;
        }
    }
}

#endif