

// #define FUNC_TRACE 1 

#ifdef FUNC_TRACE
    #define DEBUG_FUNC_START std::cout << __FILE__ << "::" <<__FUNCTION__ << " started" << std::endl;
    #define DEBUG_FUNC_END std::cout << __FILE__ << "::" <<__FUNCTION__ << " finished" << std::endl;
#else 
    #define DEBUG_FUNC_START
    #define DEBUG_FUNC_END
#endif
