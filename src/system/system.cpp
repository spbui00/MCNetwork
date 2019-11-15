#include "system.h"
#include "debug.h"


System::System(std::shared_ptr<ParameterStorage> _parameterStorage){
    DEBUG_FUNC_START

    parameterStorage=_parameterStorage;


    acceptorNumber=parameterStorage->parameters.at("acceptorNumber");


    DEBUG_FUNC_END
}

void System::setup(){
    DEBUG_FUNC_START

    //setup matrices
    pairEnergies = new double*[acceptorNumber];
    distances = new double*[acceptorNumber];
    deltaEnergies = new double*[acceptorNumber];
    for(int i=0;i<acceptorNumber;i++){
        pairEnergies[i] = new double[acceptorNumber];
        distances[i] = new double[acceptorNumber];
        deltaEnergies[i] = new double[acceptorNumber];
        pairEnergies[i][i]=0;
        distances[i][i]=0;
        deltaEnergies[i][i]=0;
    }


    // create dopants, setup position
    for(int i=0;i<acceptorNumber;i++){
        dopants.push_back(Dopant(parameterStorage->parameters.at("lenX")*enhance::random_double(0,1),parameterStorage->parameters.at("lenY")*enhance::random_double(0,1)));
    }

    // set start occupation
    std::vector<int> indicesUnoccupied {};
    for(int i=0;i<acceptorNumber;i++){indicesUnoccupied.push_back(i);}
    for(int i=0;i<acceptorNumber-parameterStorage->parameters.at("donorNumber");i++){
        int index=enhance::random_int(0,acceptorNumber-i-1);
        dopants[indicesUnoccupied[index]].occupied=true;
        indicesUnoccupied.erase(indicesUnoccupied.begin()+index);
    }

    // calc donor positions
    double donorPositions[int(parameterStorage->parameters.at("donorNumber"))][2]; 
    for(int i=0;i<parameterStorage->parameters.at("donorNumber");i++){
        donorPositions[i][0]=parameterStorage->parameters.at("lenX")*enhance::random_double(0,1);
        donorPositions[i][1]=parameterStorage->parameters.at("lenY")*enhance::random_double(0,1);
    }

    // set const energy part
    for(int i=0;i<acceptorNumber;i++){
        for(int j=0;j<parameterStorage->parameters.at("donorNumber");j++){
            dopants[i].constEnergyPart+=parameterStorage->parameters.at("I0")*1/std::sqrt(std::pow(dopants[i].posX-donorPositions[j][0],2)+std::pow(dopants[i].posY-donorPositions[j][1],2));
        }
    }

    // calc distances and pair energies
    for(int i=0;i<acceptorNumber;i++){
        for(int j=0;j<i;j++){
            pairEnergies[i][j]=parameterStorage->parameters.at("I0")/std::sqrt(std::pow(dopants[i].posX-dopants[j].posX,2)+std::pow(dopants[i].posY-dopants[j].posY,2));
            pairEnergies[j][i]=parameterStorage->parameters.at("I0")/std::sqrt(std::pow(dopants[i].posX-dopants[j].posX,2)+std::pow(dopants[i].posY-dopants[j].posY,2));
            distances[i][j]=std::sqrt(std::pow(dopants[i].posX-dopants[j].posX,2)+std::pow(dopants[i].posY-dopants[j].posY,2));
            distances[j][i]=std::sqrt(std::pow(dopants[i].posX-dopants[j].posX,2)+std::pow(dopants[i].posY-dopants[j].posY,2));
        }
    }


    // for(int i=0;i<acceptorNumber;i++){
    //     for(int j=0;j<acceptorNumber;j++){
    //         std::cout<<pairEnergies[i][j]<<"   ";
    //     }
    //     std::cout<<std::endl;
    // }
    // notOccVec = new int[acceptorNumber];
    // for(int i=0;i<acceptorNumber;i++){
    //     if (not dopants[i].occupied){
    //         notOccVec[i]=1;
    //     }
    //     else{
    //         notOccVec[i]=0;
    //     }
    // }
    

    DEBUG_FUNC_END
}

void System::calcEnergies(){
    DEBUG_FUNC_START
    // calc dopant energies
    for(int i=0;i<acceptorNumber;i++){
        dopants[i].energy=dopants[i].constEnergyPart;
    }
    for(int i=0;i<acceptorNumber;i++){
        if (not dopants[i].occupied){
            for(int j=0;j<acceptorNumber;j++){
                dopants[j].energy-=pairEnergies[i][j];
            }
        }
    }

    // Calc delta energie
    for(int i=0;i<acceptorNumber;i++){
        for(int j=0;j<acceptorNumber;j++){
            if (dopants[i].occupied && !dopants[j].occupied){
                deltaEnergies[i][j]=dopants[j].energy-dopants[i].energy-pairEnergies[i][j];
            }
            else{
                deltaEnergies[i][j]=0;
            }
        }
    }

    DEBUG_FUNC_END
}

