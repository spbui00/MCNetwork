#include "system.h"
#include "debug.h"


System::System(std::shared_ptr<ParameterStorage> parameterStorage) : parameterStorage(parameterStorage) {
    DEBUG_FUNC_START

    acceptorNumber=this->parameterStorage->parameters.at("acceptorNumber");
    hoppingSiteNumber=this->parameterStorage->parameters.at("hoppingSiteNumber");


    DEBUG_FUNC_END
}
void System::initilizeMatrices(){
    DEBUG_FUNC_START
    //setup matrices
    pairEnergies = new double*[hoppingSiteNumber];
    distances = new double*[hoppingSiteNumber];
    deltaEnergies = new double*[hoppingSiteNumber];
    for(int i=0;i<hoppingSiteNumber;i++){
        pairEnergies[i] = new double[hoppingSiteNumber];
        distances[i] = new double[hoppingSiteNumber];
        deltaEnergies[i] = new double[hoppingSiteNumber];
        pairEnergies[i][i]=0;
        distances[i][i]=0;
        deltaEnergies[i][i]=0;
    }

    donorPositions = new double*[int(parameterStorage->parameters.at("donorNumber"))];
    for(int i=0;i<parameterStorage->parameters.at("donorNumber");i++){
        donorPositions[i] = new double[2];
    }
    DEBUG_FUNC_END
}

void System::createRandomNewDevice(){
    DEBUG_FUNC_START

    // create acceptors, setup position (first part of hoppingSites, rest are electrodes)
    for(int i=0;i<acceptorNumber;i++){
        hoppingSites.push_back(new Dopant(parameterStorage->parameters.at("lenX")*enhance::random_double(0,1),parameterStorage->parameters.at("lenY")*enhance::random_double(0,1)));
    }

    // set donor positions
    for(int i=0;i<parameterStorage->parameters.at("donorNumber");i++){
        donorPositions[i][0]=parameterStorage->parameters.at("lenX")*enhance::random_double(0,1);
        donorPositions[i][1]=parameterStorage->parameters.at("lenY")*enhance::random_double(0,1);
    }

    //save device
    std::string deviceFileName=parameterStorage->workingDirecotry + "device.txt";
    ofstream deviceFile;
    deviceFile.open (deviceFileName,ios::trunc);
    deviceFile<<"acceptors: posX, posY"<<std::endl;
    for(int i=0;i<acceptorNumber;i++){
        deviceFile<<hoppingSites[i]->posX*parameterStorage->parameters["R"]<<" "<<hoppingSites[i]->posY*parameterStorage->parameters["R"]<<std::endl;
    }
    deviceFile<<std::endl;


    deviceFile<<"donors: posX, posY"<<std::endl;
    for(int i=0;i<parameterStorage->parameters.at("donorNumber");i++){
        deviceFile<<donorPositions[i][0]*parameterStorage->parameters["R"]<<" "<<donorPositions[i][1]*parameterStorage->parameters["R"]<<std::endl;
    }

    deviceFile.close();

    DEBUG_FUNC_END
}

void System::loadDevice(){
    DEBUG_FUNC_START

    std::string deviceFileName=parameterStorage->workingDirecotry + "device.txt";
    std::ifstream deviceFile(deviceFileName);
    std::string line;
    double posXBuffer, posYBuffer;

    //trahs first line
    std::getline(deviceFile, line);

    //load acceptors
    for(int i=0;i<acceptorNumber;i++){
        std::getline(deviceFile, line);
        std::istringstream iss(line);
        if(!(iss>>posXBuffer>>posYBuffer)) throw std::invalid_argument( "cant read acceptor in line: " + line);
        hoppingSites.push_back(new Dopant(posXBuffer/parameterStorage->parameters["R"],posYBuffer/parameterStorage->parameters["R"]));
    }

    //trash 2 lines
    std::getline(deviceFile, line);
    std::getline(deviceFile, line);


    // set donor positions
    for(int i=0;i<parameterStorage->parameters.at("donorNumber");i++){
        std::getline(deviceFile, line);
        std::istringstream iss(line);
        if(!(iss>>posXBuffer>>posYBuffer)) throw std::invalid_argument( "cant read donor in line: " + line);
        donorPositions[i][0]=posXBuffer/parameterStorage->parameters["R"];
        donorPositions[i][1]=posYBuffer/parameterStorage->parameters["R"];
    }
    
    deviceFile.close();

    DEBUG_FUNC_END
}

void System::getReadyForRun(){
    DEBUG_FUNC_START

    // set start occupation of acceptors
    std::vector<int> indicesUnoccupied {};
    for(int i=0;i<acceptorNumber;i++){indicesUnoccupied.push_back(i);}
    for(int i=0;i<acceptorNumber-parameterStorage->parameters.at("donorNumber");i++){
        int index=enhance::random_int(0,acceptorNumber-i-1);
        hoppingSites[indicesUnoccupied[index]]->setOccupation(true);
        indicesUnoccupied.erase(indicesUnoccupied.begin()+index);
    }




    //set electrodes
    for(ElectrodeParameters & el: parameterStorage->electrodes){
        switch (el.edge){
            case 0:
                hoppingSites.push_back(new Electrode(0,parameterStorage->parameters.at("lenY")*el.pos));
                break;
            case 1:
                hoppingSites.push_back(new Electrode(parameterStorage->parameters.at("lenX"),parameterStorage->parameters.at("lenY")*el.pos));
                break;
            case 2:
                hoppingSites.push_back(new Electrode(parameterStorage->parameters.at("lenX")*el.pos,0));
                break;
            case 3:
                hoppingSites.push_back(new Electrode(parameterStorage->parameters.at("lenX")*el.pos,parameterStorage->parameters.at("lenY")));
                break;
        }
    }





    //solve laplace eq
    finEle= new FiniteElemente(parameterStorage->parameters.at("lenX"),parameterStorage->parameters.at("lenY"),parameterStorage->parameters.at("finiteElementsResolution"));
    //set electrodes
    for(ElectrodeParameters & el: parameterStorage->electrodes){
        finEle->setElectrode(parameterStorage->parameters.at("lenX")*el.pos-0.5*parameterStorage->parameters.at("electrodeWidth"),parameterStorage->parameters.at("lenX")*el.pos+0.5*parameterStorage->parameters.at("electrodeWidth"),el.edge,el.voltage);
    }
    finEle->initRun();
    finEle->run();



    // calc distances and pair energies
    for(int i=0;i<hoppingSiteNumber;i++){
        for(int j=0;j<i;j++){
            pairEnergies[i][j]=parameterStorage->parameters.at("I0")/std::sqrt(std::pow(hoppingSites[i]->posX-hoppingSites[j]->posX,2)+std::pow(hoppingSites[i]->posY-hoppingSites[j]->posY,2));
            pairEnergies[j][i]=parameterStorage->parameters.at("I0")/std::sqrt(std::pow(hoppingSites[i]->posX-hoppingSites[j]->posX,2)+std::pow(hoppingSites[i]->posY-hoppingSites[j]->posY,2));
            distances[i][j]=std::sqrt(std::pow(hoppingSites[i]->posX-hoppingSites[j]->posX,2)+std::pow(hoppingSites[i]->posY-hoppingSites[j]->posY,2));
            distances[j][i]=std::sqrt(std::pow(hoppingSites[i]->posX-hoppingSites[j]->posX,2)+std::pow(hoppingSites[i]->posY-hoppingSites[j]->posY,2));
            // std::cout<<"pair e "<<pairEnergies[i][j]<<" I0 "<<parameterStorage->parameters.at("I0")<<std::endl;
        }
    }

    // set const energy part
    //donor interaction
    for(int i=0;i<acceptorNumber;i++){
        for(int j=0;j<parameterStorage->parameters.at("donorNumber");j++){
            hoppingSites[i]->constEnergyPart+=parameterStorage->parameters.at("I0")*1/std::sqrt(std::pow(hoppingSites[i]->posX-donorPositions[j][0],2)+std::pow(hoppingSites[i]->posY-donorPositions[j][1],2));
        }
    }
    //potential
    for(int i=0;i<hoppingSiteNumber;i++){
        // std::cout<<i<<" donor interaction "<< hoppingSites[i]->constEnergyPart <<" pot "<< finEle->getPotential(hoppingSites[i]->posX,hoppingSites[i]->posY)*parameterStorage->parameters.at("e")/parameterStorage->parameters.at("kT")<<std::endl;
        hoppingSites[i]->constEnergyPart+=finEle->getPotential(hoppingSites[i]->posX,hoppingSites[i]->posY)*parameterStorage->parameters.at("e")/parameterStorage->parameters.at("kT");
    }

    //set electrode energy
    for(int i=acceptorNumber;i<hoppingSiteNumber;i++){
        hoppingSites[i]->energy=hoppingSites[i]->constEnergyPart;
    }

    // -------- new -------------
    for(int i=0;i<acceptorNumber;i++){
        hoppingSites[i]->energy2=hoppingSites[i]->constEnergyPart;
    }

    DEBUG_FUNC_END
}



void System::calcEnergies(){
    DEBUG_FUNC_START

    // -------- new -------------

    for(int i=0;i<acceptorNumber;i++){
        hoppingSites[i]->energy2=hoppingSites[i]->constEnergyPart;
    }







    // -------- old -------------
    // calc dopant energies
    for(int i=0;i<acceptorNumber;i++){
        hoppingSites[i]->energy=hoppingSites[i]->constEnergyPart;
    }
    for(int i=0;i<acceptorNumber;i++){ //coulomb interaction only with acceptors and..
        if (not hoppingSites[i]->getOccupation()){ //.. only if they are unoccupied
            for(int j=0;j<acceptorNumber;j++){
                hoppingSites[j]->energy-=pairEnergies[i][j];
            }
        }
    }

    // Calc delta energie
    //acceptor acceptor hopp
    for(int i=0;i<acceptorNumber;i++){
        if(hoppingSites[i]->getOccupation()){
            for(int j=0;j<acceptorNumber;j++){ 
                if (!hoppingSites[j]->getOccupation()){
                    deltaEnergies[i][j]=hoppingSites[j]->energy-hoppingSites[i]->energy-pairEnergies[i][j];
                }
                else{
                    deltaEnergies[i][j]=0;
                }
            }
        }
        else{
            for(int j=0;j<acceptorNumber;j++){ 
                deltaEnergies[i][j]=0;
            }
        }
    }

    //acceptor electrode hopp
    for(int i=0;i<acceptorNumber;i++){
        if (hoppingSites[i]->getOccupation()){
            for(int j=acceptorNumber;j<hoppingSiteNumber;j++){ 
                deltaEnergies[i][j]=hoppingSites[j]->energy-hoppingSites[i]->energy;
            }
        }
        else{
            for(int j=acceptorNumber;j<hoppingSiteNumber;j++){ 
                deltaEnergies[i][j]=0;
            }
        }
    }
    //electrode acceptor hopp
    for(int j=0;j<acceptorNumber;j++){ 
        if (!hoppingSites[j]->getOccupation()){
            for(int i=acceptorNumber;i<hoppingSiteNumber;i++){
                deltaEnergies[i][j]=hoppingSites[j]->energy-hoppingSites[i]->energy;
            }
        }
        else{
            for(int i=acceptorNumber;i<hoppingSiteNumber;i++){
                deltaEnergies[i][j]=0;
            }
        }
    }
    //electrode electrode hopp
    for(int i=acceptorNumber;i<hoppingSiteNumber;i++){
        for(int j=acceptorNumber;j<i;j++){ 
            deltaEnergies[i][j]=hoppingSites[j]->energy-hoppingSites[i]->energy;
            deltaEnergies[j][i]=hoppingSites[i]->energy-hoppingSites[j]->energy;
        }
    }

    // for(int i=0;i<hoppingSiteNumber;i++){
    //     for(int j=0;j<hoppingSiteNumber;j++){
    //         std::cout<<deltaEnergies[i][j]<<" ";
    //     }
    //     std::cout<<std::endl;
    // }
    // std::cout<<std::endl;



    DEBUG_FUNC_END
}

void System::updatePotential(){
    DEBUG_FUNC_START
    
    //reset old potential
    for(int i=0;i<hoppingSiteNumber;i++){
        hoppingSites[i]->constEnergyPart-=finEle->getPotential(hoppingSites[i]->posX,hoppingSites[i]->posY)*parameterStorage->parameters.at("e")/parameterStorage->parameters.at("kT");
    }

    //recalc potential
    for(int i=0;i < parameterStorage->electrodes.size();i++){
        finEle->updateElectrodeVoltage(i,parameterStorage->electrodes[i].voltage);
    }
    finEle->run();

    //set new potential
    for(int i=0;i<hoppingSiteNumber;i++){
        hoppingSites[i]->constEnergyPart+=finEle->getPotential(hoppingSites[i]->posX,hoppingSites[i]->posY)*parameterStorage->parameters.at("e")/parameterStorage->parameters.at("kT");
    }
    DEBUG_FUNC_END
}



void System::increaseTime(double const & ratesSum){
    DEBUG_FUNC_START

    time+=std::log(enhance::random_double(0,1))/(-1*ratesSum);

    DEBUG_FUNC_END
}

std::string System::getState(){
    DEBUG_FUNC_START

    std::string occupationState="";
    for(int i=0;i<acceptorNumber;i++){
        occupationState+=hoppingSites[i]->getOccupation() ? "1" : "0";
    }

    // std::cout<<occupationState<<std::endl;
    return occupationState;

    DEBUG_FUNC_END
}
