#pragma once

#ifndef NDEBUG
#include "enhance.hpp"
#include <iostream>
#include <set>
#include <string>

class myDebugger {
public:
    static int funcTraceLevel;
    static void debugPrint(std::string message, std::string key);
    static std::set<std::string> debugKeys;
};

#define DEBUG_FUNC_START                                                                                                                          \
    if (myDebugger::debugKeys.find("funcTrace") != myDebugger::debugKeys.end()) {                                                                 \
        myDebugger::funcTraceLevel++;                                                                                                             \
        std::cout << enhance::multipliplyString("  ", myDebugger::funcTraceLevel) << __FILE__ << "::" << __FUNCTION__ << " started" << std::endl; \
    }
#define DEBUG_FUNC_END                                                                                                                             \
    if (myDebugger::debugKeys.find("funcTrace") != myDebugger::debugKeys.end()) {                                                                  \
        std::cout << enhance::multipliplyString("  ", myDebugger::funcTraceLevel) << __FILE__ << "::" << __FUNCTION__ << " finished" << std::endl; \
        myDebugger::funcTraceLevel--;                                                                                                              \
    }
#define DBEUG_LOG(message, key) myDebugger::debugPrint(message, key);
#else
#define DEBUG_FUNC_START
#define DEBUG_FUNC_END
#define DBEUG_LOG(message, key)
#endif
