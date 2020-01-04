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
        if(!(iss>>name)) throw std::invalid_argument( "cant read line: " + line);
        if(name == "electrode"){
            double pos, voltage;
            int edge;
            if(!(iss>>pos>>edge>>voltage)) throw std::invalid_argument( "cant read electrode: " + line);
            electrodes.push_back({pos,edge,voltage});
        }
        else if (name == "gate")
        {
            if(!(iss>>gate)) throw std::invalid_argument( "cant read gate: " + line);
        }
        else{
            if(!(iss>>val)) throw std::invalid_argument( "can't read line: " + line);
            parameters[name]=val;    
            // std::cout<<name<<parameters[name]<<std::endl;
        }
    }


    parameters["kT"]=parameters["k"]*parameters["T"];
    parameters["donorNumber"]=int(parameters["acceptorNumber"]*parameters["compensationFactor"]);
    parameters["hoppingSiteNumber"]=parameters["acceptorNumber"]+electrodes.size();

    // convert lens in dimensions of R    
    parameters["R"]=std::pow(parameters["lenX"]*parameters["lenY"]/parameters["acceptorNumber"],0.5);
    parameters["lenX"]=parameters["lenX"]/parameters["R"];
    parameters["lenY"]=parameters["lenY"]/parameters["R"];
    parameters["a"]=parameters["a"]/parameters["R"];
    parameters["electrodeWidth"]=parameters["electrodeWidth"]/parameters["R"];

    std::cout<<"R = "<<parameters["R"]<<"nm, a/R = "<<parameters["a"]<<" I0 = "<<parameters["I0"]<<" meV == "<<parameters["I0"]*0.001*parameters["e"]/parameters["kT"]<<"kT corresponds to eps = "<< parameters["e"]*parameters["e"]/(4*M_PI*parameters["eps0"]*parameters["R"]*1e-9*parameters["I0"]*0.001*parameters["e"]) <<std::endl;

    // convert energies in dimensions of kT
    parameters["I0"]=parameters["I0"]*0.001*parameters["e"]/parameters["kT"]; //I0 given in meV


  

    DEBUG_FUNC_END
}
