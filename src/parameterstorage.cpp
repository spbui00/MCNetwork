#include "parameterstorage.h"
ParameterStorage::ParameterStorage(std::string filename)
{
    DEBUG_FUNC_START

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


    parameters["kT"]=parameters["k"]*parameters["T"];
    parameters["donorNumber"]=int(parameters["acceptorNumber"]*parameters["compensationFactor"]);
    
    // convert lens in dimensions of R    
    parameters["R"]=std::pow(parameters["lenX"]*parameters["lenY"]/parameters["acceptorNumber"],0.5);
    parameters["lenX"]=parameters["lenX"]/parameters["R"];
    parameters["lenY"]=parameters["lenY"]/parameters["R"];
    parameters["a"]=parameters["a"]/parameters["R"];
    // convert energies in dimensions of kT
    parameters["I0"]=parameters["I0"]*0.001*parameters["e"]/parameters["kT"];



  

    DEBUG_FUNC_END
}
