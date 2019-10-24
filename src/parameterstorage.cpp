#include "parameterstorage.h"
#ifndef NDEBUG
   #include "debug.h"
#endif
ParameterStorage::ParameterStorage(std::string filename)
{
    #ifndef NDEBUG
        DEBUG_FUNC_START
    #endif
    std::ifstream infile(filename);
    std::string line;
    double val;  //buffer for values
    std::string name; //buffer for strings

    while (std::getline(infile, line))
    {   
            // std::cout<<line<<std::endl;

        //remove comments
        auto pos=line.find("#");
        if (pos==0 or line.size()==0) continue;
        if (pos==std::string::npos);
        else
        {
            line.erase(pos,std::string::npos);
        }

        std::istringstream iss(line);
        if(!(iss>>name)) throw std::invalid_argument( "cant read line: "+line );
        {
            if(!(iss>>val)) throw std::invalid_argument( "can't read line: "+line );
            parameters[name]=val;    
            // std::cout<<name<<parameters[name]<<std::endl;
        }
    }
  

    #ifndef NDEBUG
        DEBUG_FUNC_END
    #endif    
}
