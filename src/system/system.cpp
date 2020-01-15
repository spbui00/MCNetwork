#include "system.h"
#include "debug.h"


System::System(std::shared_ptr<ParameterStorage> parameterStorage) : parameterStorage(parameterStorage) {
    DEBUG_FUNC_START

    acceptorNumber    = parameterStorage->parameters.at("acceptorNumber");
    hoppingSiteNumber = parameterStorage->parameters.at("hoppingSiteNumber");
    locLenA           = parameterStorage->parameters.at("a");
    electrodeNumber   = hoppingSiteNumber-acceptorNumber;


    DEBUG_FUNC_END
}
void System::initilizeMatrices(){
    DEBUG_FUNC_START
    //setup matrices
    pairEnergies  = new double[acceptorNumber*acceptorNumber];
    distances     = new double[hoppingSiteNumber*hoppingSiteNumber];
    deltaEnergies = new double[hoppingSiteNumber*hoppingSiteNumber];
    energies       = new double[hoppingSiteNumber];
    currentCounter = new double[hoppingSiteNumber];
    rates = std::make_shared<std::vector<double>>(hoppingSiteNumber*hoppingSiteNumber);

    donorPositionsX     = new double[int(parameterStorage->parameters.at("donorNumber"))];
    donorPositionsY     = new double[int(parameterStorage->parameters.at("donorNumber"))];
    acceptorPositionsX  = new double[acceptorNumber];
    acceptorPositionsY  = new double[acceptorNumber];
    electrodePositionsX = new double[electrodeNumber];
    electrodePositionsY = new double[electrodeNumber];

    occupation = new bool[acceptorNumber];
    for(int i=0;i<acceptorNumber;i++){
        occupation[i]=false;
        // currentState += "0";
    }


    DEBUG_FUNC_END
}

void System::createRandomNewDevice(){
    DEBUG_FUNC_START

    // create acceptors, setup position (first part of hoppingSites, rest are electrodes)
    for(int i=0;i<acceptorNumber;i++){
        acceptorPositionsX[i]=parameterStorage->parameters.at("lenX")*enhance::random_double(0,1);
        acceptorPositionsY[i]=parameterStorage->parameters.at("lenY")*enhance::random_double(0,1);
    }

    // set donor positions
    for(int i=0;i<parameterStorage->parameters.at("donorNumber");i++){
        donorPositionsX[i]=parameterStorage->parameters.at("lenX")*enhance::random_double(0,1);
        donorPositionsY[i]=parameterStorage->parameters.at("lenY")*enhance::random_double(0,1);
    }

    //save device
    std::string deviceFileName=parameterStorage->workingDirecotry + "device.txt";
    ofstream deviceFile;
    deviceFile.open (deviceFileName,ios::trunc);
    deviceFile<<"acceptors: posX, posY"<<std::endl;
    for(int i=0;i<acceptorNumber;i++){
        deviceFile<<acceptorPositionsX[i]*parameterStorage->parameters["R"]<<" "<<acceptorPositionsY[i]*parameterStorage->parameters["R"]<<std::endl;
    }
    deviceFile<<std::endl;


    deviceFile<<"donors: posX, posY"<<std::endl;
    for(int i=0;i<parameterStorage->parameters.at("donorNumber");i++){
        deviceFile<<donorPositionsX[i]*parameterStorage->parameters["R"]<<" "<<donorPositionsY[i]*parameterStorage->parameters["R"]<<std::endl;
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

    //trash first line
    std::getline(deviceFile, line);

    //load acceptors
    for(int i=0;i<acceptorNumber;i++){
        std::getline(deviceFile, line);
        std::istringstream iss(line);
        if(!(iss>>posXBuffer>>posYBuffer)) throw std::invalid_argument( "cant read acceptor in line: " + line);
        acceptorPositionsX[i]=posXBuffer/parameterStorage->parameters["R"];
        acceptorPositionsY[i]=posYBuffer/parameterStorage->parameters["R"];
    }

    //trash 2 lines
    std::getline(deviceFile, line);
    std::getline(deviceFile, line);


    // set donor positions
    for(int i=0;i<parameterStorage->parameters.at("donorNumber");i++){
        std::getline(deviceFile, line);
        std::istringstream iss(line);
        if(!(iss>>posXBuffer>>posYBuffer)) throw std::invalid_argument( "cant read donor in line: " + line);
        donorPositionsX[i]=posXBuffer/parameterStorage->parameters["R"];
        donorPositionsY[i]=posYBuffer/parameterStorage->parameters["R"];
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
        occupation[indicesUnoccupied[index]]=true;
        hasedCurrentState+=enhance::fastExp2(indicesUnoccupied[index]); //compute hash
        indicesUnoccupied.erase(indicesUnoccupied.begin()+index);
    }


    //set electrodes
    for(int i=0; i< electrodeNumber; i++){
        switch (parameterStorage->electrodes[i].edge){
            case 0:
                electrodePositionsX[i]=0;
                electrodePositionsY[i]=parameterStorage->parameters.at("lenY")*parameterStorage->electrodes[i].pos;
                break;
            case 1:
                electrodePositionsX[i]=parameterStorage->parameters.at("lenX");
                electrodePositionsY[i]=parameterStorage->parameters.at("lenY")*parameterStorage->electrodes[i].pos;
                break;
            case 2:
                electrodePositionsX[i]=parameterStorage->parameters.at("lenX")*parameterStorage->electrodes[i].pos;
                electrodePositionsY[i]=0;
                break;
            case 3:
                electrodePositionsX[i]=parameterStorage->parameters.at("lenX")*parameterStorage->electrodes[i].pos;
                electrodePositionsY[i]=parameterStorage->parameters.at("lenY");
                break;
        }
    }





    //solve laplace eq
    finEle= new FiniteElemente(parameterStorage->parameters.at("lenX"),parameterStorage->parameters.at("lenY"),parameterStorage->parameters.at("finiteElementsResolution"));
    //set electrodes
    for(int i=0; i< electrodeNumber; i++){
        switch (parameterStorage->electrodes[i].edge){
            case 0:
            case 1:
                finEle->setElectrode(parameterStorage->parameters.at("lenY")*parameterStorage->electrodes[i].pos-0.5*parameterStorage->parameters.at("electrodeWidth"),parameterStorage->parameters.at("lenY")*parameterStorage->electrodes[i].pos+0.5*parameterStorage->parameters.at("electrodeWidth"),parameterStorage->electrodes[i].edge,parameterStorage->electrodes[i].voltage);
                break;
            case 2:
            case 3:
                finEle->setElectrode(  parameterStorage->parameters.at("lenX")*parameterStorage->electrodes[i].pos-0.5*parameterStorage->parameters.at("electrodeWidth"),parameterStorage->parameters.at("lenX")*parameterStorage->electrodes[i].pos+0.5*parameterStorage->parameters.at("electrodeWidth"),parameterStorage->electrodes[i].edge,parameterStorage->electrodes[i].voltage);
                break;
        }
    }

    finEle->initRun();
    finEle->run();



    // calc distances and pair energies
    // acc<->acc
    for(int i=0;i<acceptorNumber;i++){
        for(int j=0;j<i;j++){
            distances[i*hoppingSiteNumber+j]=std::sqrt(std::pow(acceptorPositionsX[i]-acceptorPositionsX[j],2)+std::pow(acceptorPositionsY[i]-acceptorPositionsY[j],2));
            distances[j*hoppingSiteNumber+i]=std::sqrt(std::pow(acceptorPositionsX[i]-acceptorPositionsX[j],2)+std::pow(acceptorPositionsY[i]-acceptorPositionsY[j],2));
            pairEnergies[i*acceptorNumber+j]=-parameterStorage->parameters.at("I0")/distances[i*hoppingSiteNumber+j];
            pairEnergies[j*acceptorNumber+i]=-parameterStorage->parameters.at("I0")/distances[j*hoppingSiteNumber+i];
            // std::cout<<"pair e "<<pairEnergies[i*hoppingSiteNumber+j]<<std::endl;
        }
        pairEnergies[i*(acceptorNumber+1)] = 0; //[i*(hoppingSiteNumber+1)] = [i,i]
        distances   [i*(hoppingSiteNumber+1)] = 0;
    }
    // acc->el
    for(int i=0;i<acceptorNumber;i++){
        for(int j=0;j<electrodeNumber;j++){
            distances[i*hoppingSiteNumber+(j+acceptorNumber)]=std::sqrt(std::pow(acceptorPositionsX[i]-electrodePositionsX[j],2)+std::pow(acceptorPositionsY[i]-electrodePositionsY[j],2));
        }
    }
    // el->acc
    for(int i=0;i<electrodeNumber;i++){
        for(int j=0;j<acceptorNumber;j++){
            distances[(i+acceptorNumber)*hoppingSiteNumber+j]=std::sqrt(std::pow(electrodePositionsX[i]-acceptorPositionsX[j],2)+std::pow(electrodePositionsY[i]-acceptorPositionsY[j],2));
        }
    }
    // el<->el
    for(int i=0;i<electrodeNumber;i++){
        for(int j=0;j<i;j++){
            distances[(i+acceptorNumber)*hoppingSiteNumber+(j+acceptorNumber)]=std::sqrt(std::pow(electrodePositionsX[i]-electrodePositionsX[j],2)+std::pow(electrodePositionsY[i]-electrodePositionsY[j],2));
            distances[(j+acceptorNumber)*hoppingSiteNumber+(i+acceptorNumber)]=std::sqrt(std::pow(electrodePositionsX[i]-electrodePositionsX[j],2)+std::pow(electrodePositionsY[i]-electrodePositionsY[j],2));
        }
        distances[(i+acceptorNumber)*(hoppingSiteNumber+1)]=0;
    }


    // set acceptor energies
    // donor interaction
    for(int i=0;i<acceptorNumber;i++){
        energies[i]=0;
        for(int j=0;j<parameterStorage->parameters.at("donorNumber");j++){
            energies[i]+=parameterStorage->parameters.at("I0")*1/std::sqrt(std::pow(acceptorPositionsX[i]-donorPositionsX[j],2)+std::pow(acceptorPositionsY[i]-donorPositionsY[j],2));
        }
    }
    //potential
    for(int i=0;i<acceptorNumber;i++){
        energies[i]+=finEle->getPotential(acceptorPositionsX[i],acceptorPositionsY[i])*parameterStorage->parameters.at("e")/parameterStorage->parameters.at("kT");
        // std::cout<<i<<" donor interaction "<< hoppingSites[i]->constEnergyPart <<" pot "<< finEle->getPotential(hoppingSites[i]->posX,hoppingSites[i]->posY)*parameterStorage->parameters.at("e")/parameterStorage->parameters.at("kT")<<std::endl;
    }
    //set coulomb part (with start occupation)
    for(int i=0;i<acceptorNumber;i++){ //coulomb interaction only with acceptors and..
        if (not occupation[i]){ //.. only if they are unoccupied
            for(int j=0;j<acceptorNumber;j++){ //self interaction is ignored, bc pairEnergies[i][i]=0
                energies[j]+=pairEnergies[i*acceptorNumber+j];
            }
        }
    }

    //set electrode energy
    for(int i=0;i<electrodeNumber;i++){
        energies[i+acceptorNumber]=finEle->getPotential(electrodePositionsX[i],electrodePositionsY[i])*parameterStorage->parameters.at("e")/parameterStorage->parameters.at("kT");
    }


    //set deltaEnergies and rates (only constant part = el-el interaction)
    for(int i=acceptorNumber;i<hoppingSiteNumber;i++){
        (*rates)[i*(acceptorNumber+1)]=0;
        for(int j=acceptorNumber;j<i;j++){
            deltaEnergies[i*hoppingSiteNumber+j]=energies[j]-energies[i];
            deltaEnergies[j*hoppingSiteNumber+i]=energies[i]-energies[j];
            // std::cout<<"deltaEnergies el el "<< distances[i*hoppingSiteNumber+j]<<std::endl;

            if (deltaEnergies[i*hoppingSiteNumber+j] < 0){
                (*rates)[i*hoppingSiteNumber+j]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA);
                (*rates)[j*hoppingSiteNumber+i]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA-deltaEnergies[j*hoppingSiteNumber+i]);
            }
            else if (deltaEnergies[i*hoppingSiteNumber+j] > 0){
                (*rates)[i*hoppingSiteNumber+j]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA-deltaEnergies[i*hoppingSiteNumber+j]);
                (*rates)[j*hoppingSiteNumber+i]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA);
            }
            else{
                (*rates)[i*hoppingSiteNumber+j]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA);
                (*rates)[j*hoppingSiteNumber+i]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA);
            }
        }
    }
    DEBUG_FUNC_END
}

void System::updateRates(bool storeKnowStates /* = false*/){
    DEBUG_FUNC_START


    if (knownRates.count(hasedCurrentState)){
        rates    = knownRates   .at(hasedCurrentState);
        ratesSum = knownRatesSum.at(hasedCurrentState);
        // std::cout<<"found state: "<<state<<std::endl;
        
    }
    else{
        ratesSum=0;


        if (rates.use_count() > 1){
            rates = std::make_shared<std::vector<double>>(hoppingSiteNumber*hoppingSiteNumber);
        }

        //acc acc hopp
        for(int i=0;i<acceptorNumber;i++){
            if (occupation[i]){
                for(int j=0;j<acceptorNumber;j++){
                    if (not occupation[j]){
                        // std::cout<<" ----- "<<i<<" "<<j<<std::endl;

                        deltaEnergies[i*hoppingSiteNumber+j]=energies[j]-energies[i]+pairEnergies[i*acceptorNumber+j];

                        // std::cout<<" ei "<<energies[i]<<" ej "<<energies[j]<<" epair "<<pairEnergies[i*acceptorNumber+j]<<std::endl;

                        if (deltaEnergies[i*hoppingSiteNumber+j] <= 0){
                            (*rates)[i*hoppingSiteNumber+j]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA);
                        }
                        else{
                            (*rates)[i*hoppingSiteNumber+j]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA-deltaEnergies[i*hoppingSiteNumber+j]);
                        }
                        ratesSum+=(*rates)[i*hoppingSiteNumber+j];
                    }
                    else{
                        (*rates)[i*hoppingSiteNumber+j]=0;
                        // deltaEnergies[i*hoppingSiteNumber+j]=0;
                    }
                }
            }
            else{
                for(int j=0;j<acceptorNumber;j++){ 
                    (*rates)[i*hoppingSiteNumber+j]=0;
                    // deltaEnergies[i*hoppingSiteNumber+j]=0;
                }
            }
        }

        // el-acc hopp
        for(int i=acceptorNumber;i<hoppingSiteNumber;i++){
            for(int j=0;j<acceptorNumber;j++){
                if (not occupation[j]){
                    deltaEnergies[i*hoppingSiteNumber+j]=energies[j]-energies[i];
                    if (deltaEnergies[i*hoppingSiteNumber+j] <= 0){
                        (*rates)[i*hoppingSiteNumber+j]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA);
                    }
                    else{
                        (*rates)[i*hoppingSiteNumber+j]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA-deltaEnergies[i*hoppingSiteNumber+j]);
                    }
                    ratesSum+=(*rates)[i*hoppingSiteNumber+j];
                }
                else{
                    (*rates)[i*hoppingSiteNumber+j]=0;
                    // deltaEnergies[i*hoppingSiteNumber+j]=0;
                }
            }
        }


        //acc-el hopp
        for(int i=0;i<acceptorNumber;i++){
            if (occupation[i]){
                for(int j=acceptorNumber;j<hoppingSiteNumber;j++){
                    deltaEnergies[i*hoppingSiteNumber+j]=energies[j]-energies[i];
                    if (deltaEnergies[i*hoppingSiteNumber+j] <= 0){
                        (*rates)[i*hoppingSiteNumber+j]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA);
                    }
                    else{
                        (*rates)[i*hoppingSiteNumber+j]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA-deltaEnergies[i*hoppingSiteNumber+j]);
                    }
                    ratesSum+=(*rates)[i*hoppingSiteNumber+j];
                }
            }
            else{
                for(int j=acceptorNumber;j<hoppingSiteNumber;j++){ 
                    (*rates)[i*hoppingSiteNumber+j]=0;
                    // deltaEnergies[i*hoppingSiteNumber+j]=0;
                }
            }
        }

        if (storeKnowStates){
            knownRates   [hasedCurrentState]=rates;
            knownRatesSum[hasedCurrentState]=ratesSum;
        }
    }

    DEBUG_FUNC_END
}

void System::updatePotential(){
    DEBUG_FUNC_START
    
    //reset old potential
    for(int i=0;i<acceptorNumber;i++){
        energies[i]-=finEle->getPotential(acceptorPositionsX[i],acceptorPositionsY[i])*parameterStorage->parameters.at("e")/parameterStorage->parameters.at("kT");
    }
    for(int i=0;i<electrodeNumber;i++){
        energies[(i+acceptorNumber)]-=finEle->getPotential(electrodePositionsX[i],electrodePositionsY[i])*parameterStorage->parameters.at("e")/parameterStorage->parameters.at("kT");
    }


    //recalc potential
    for(int i=0;i < parameterStorage->electrodes.size();i++){
        finEle->updateElectrodeVoltage(i,parameterStorage->electrodes[i].voltage);
    }
    finEle->run();

    //set new potential
    for(int i=0;i<acceptorNumber;i++){
        energies[i]+=finEle->getPotential(acceptorPositionsX[i],acceptorPositionsY[i])*parameterStorage->parameters.at("e")/parameterStorage->parameters.at("kT");
    }
    for(int i=0;i<electrodeNumber;i++){
        energies[(i+acceptorNumber)]+=finEle->getPotential(electrodePositionsX[i],electrodePositionsY[i])*parameterStorage->parameters.at("e")/parameterStorage->parameters.at("kT");
    }
    DEBUG_FUNC_END
}

void System::makeSwap(){
    DEBUG_FUNC_START
    double rndNumber=enhance::random_double(0,ratesSum);
    double partRatesSum=0;
    for(int i=0; i<hoppingSiteNumber;i++){
        for(int j=0; j<hoppingSiteNumber;j++){
            partRatesSum+=(*rates)[i*hoppingSiteNumber+j];
            // std::cout<<partRatesSum<<" "<<(*rates)[i*hoppingSiteNumber+j]<<" "<<rndNumber<<" "<<ratesSum<<std::endl;
            if(partRatesSum > rndNumber){
                lastSwapped1=i;
                lastSwapped2=j;

                currentCounter[i]--;
                currentCounter[j]++;

                // std::cout<<"swapped "<<i<<" "<<j<<" "<<setw(9);
                goto endDoubleLoop;
            }
        }
        // if(i== hoppingSiteNumber-1){
        //     std::cout<<"no swapp found!"<<partRatesSum<<" "<<rndNumber<<" "<<ratesSum<<" ";
        // }
    }    
    endDoubleLoop:;

    updateAfterSwap();


    // for(int i=0; i<hoppingSiteNumber;i++){
    //     // std::cout<<"i "<<i<<" curr "<<currentCounter[i]<<" "<<std::endl;;
    // }


    DEBUG_FUNC_END
}

void System::updateAfterSwap(){
    DEBUG_FUNC_START

    if (lastSwapped1 < acceptorNumber){ //last swapped1 = acceptor
        occupation[lastSwapped1]=false;
        //update energy
        for(int j=0;j<acceptorNumber;j++){
            energies[j]+=pairEnergies[lastSwapped1*acceptorNumber+j];
        }
        //update hashedState
        hasedCurrentState-=enhance::fastExp2(lastSwapped1); //compute hash
    }

    if (lastSwapped2 < acceptorNumber){ //last swapped2 = acceptor
        occupation[lastSwapped2]=true;
        for(int j=0;j<acceptorNumber;j++){
            energies[j]-=pairEnergies[lastSwapped2*acceptorNumber+j];
        }
        hasedCurrentState+=enhance::fastExp2(lastSwapped2); //compute hash
    }

    // std::cout<<hasedCurrentState<<std::endl;

    DEBUG_FUNC_END
}



void System::increaseTime(){
    DEBUG_FUNC_START

    time+=std::log(enhance::random_double(0,1))/(-1*ratesSum);

    DEBUG_FUNC_END
}
