#include "system.h"
#include "debug.h"


System::System(const std::shared_ptr<ParameterStorage> & parameterStorage) : parameterStorage(parameterStorage) {
    DEBUG_FUNC_START

    acceptorNumber    = parameterStorage->parameters.at("acceptorNumber"   );
    hoppingSiteNumber = parameterStorage->parameters.at("hoppingSiteNumber");
    locLenA           = parameterStorage->parameters.at("a"                );
    storingMode       = parameterStorage->parameters.at("storingMode"      );
    electrodeNumber   = hoppingSiteNumber-acceptorNumber;

    konwnPartRatesSumList.reset(new std::unordered_map<unsigned long long,std::shared_ptr<std::vector<double>>>());
    knownRatesSum        .reset(new std::unordered_map<unsigned long long,double>());

    storeKnownStates = new bool(true);

    mutex=std::make_shared<std::shared_mutex>();

    DEBUG_FUNC_END
}
System::System(const System & oldSys, bool shareMemory /*= true*/) : 
                                        parameterStorage     (oldSys.parameterStorage     ),
                                        acceptorNumber       (oldSys.acceptorNumber       ),
                                        hoppingSiteNumber    (oldSys.hoppingSiteNumber    ),
                                        electrodeNumber      (oldSys.electrodeNumber      ),
                                        donorPositionsX      (oldSys.donorPositionsX      ),
                                        donorPositionsY      (oldSys.donorPositionsY      ),
                                        acceptorPositionsX   (oldSys.acceptorPositionsX   ),
                                        acceptorPositionsY   (oldSys.acceptorPositionsY   ),
                                        electrodePositionsX  (oldSys.electrodePositionsX  ),
                                        electrodePositionsY  (oldSys.electrodePositionsY  ),
                                        distances            (oldSys.distances            ),
                                        pairEnergies         (oldSys.pairEnergies         ),
                                        hasedCurrentState    (oldSys.hasedCurrentState    ),
                                        locLenA              (oldSys.locLenA              ),
                                        constantRatesSumPart (oldSys.constantRatesSumPart ),
                                        finEle               (oldSys.finEle               ),
                                        storeKnownStates     (oldSys.storeKnownStates     ),
                                        konwnPartRatesSumList(oldSys.konwnPartRatesSumList),
                                        storingMode          (oldSys.storingMode          ),
                                        knownRatesSum        (oldSys.knownRatesSum        ),
                                        mutex                (oldSys.mutex                ) {
    DEBUG_FUNC_START

    if(not oldSys.readyForRun){
        throw std::invalid_argument("copying system that isnt ready for run is forbidden!");
    }

    energies       = new double[hoppingSiteNumber];
    currentCounter = new double[hoppingSiteNumber];
    for(int i=0;i<hoppingSiteNumber;i++){
        energies      [i] = oldSys.energies      [i];
        currentCounter[i] = oldSys.currentCounter[i];
    }

    occupation = new bool[acceptorNumber];
    for(int i=0;i<acceptorNumber;i++){
        occupation[i] = oldSys.occupation[i];
    }

    deltaEnergies = new double[hoppingSiteNumber*hoppingSiteNumber];
    rates         = new double[hoppingSiteNumber*hoppingSiteNumber];
    for(int i=0;i<hoppingSiteNumber*hoppingSiteNumber;i++){
        deltaEnergies[i] = oldSys.deltaEnergies[i];
        rates        [i] = oldSys.rates        [i];
    }

    if (not shareMemory){
        konwnPartRatesSumList.reset(new std::unordered_map<unsigned long long,std::shared_ptr<std::vector<double>>>());
        knownRatesSum        .reset(new std::unordered_map<unsigned long long,double>());
        // distances      = new double[hoppingSiteNumber*hoppingSiteNumber];
        // for(int i=0;i<hoppingSiteNumber*hoppingSiteNumber;i++){
        //     distances[i]=oldSys.distances[i];
        // }
        // pairEnergies   = new double[acceptorNumber*acceptorNumber];
        // for(int i=0;i<acceptorNumber*acceptorNumber;i++){
        //     pairEnergies[i]=oldSys.pairEnergies[i];
        // }
    }

    readyForRun=true;
    DEBUG_FUNC_END
}

void System::initilizeMatrices(){
    DEBUG_FUNC_START
    //setup matrices
    pairEnergies   = new double[acceptorNumber*acceptorNumber];
    distances      = new double[hoppingSiteNumber*hoppingSiteNumber];
    deltaEnergies  = new double[hoppingSiteNumber*hoppingSiteNumber];
    rates          = new double[hoppingSiteNumber*hoppingSiteNumber];
    energies       = new double[hoppingSiteNumber];
    currentCounter = new double[hoppingSiteNumber];

    donorPositionsX     = new double[int(parameterStorage->parameters.at("donorNumber"))];
    donorPositionsY     = new double[int(parameterStorage->parameters.at("donorNumber"))];
    acceptorPositionsX  = new double[acceptorNumber ];
    acceptorPositionsY  = new double[acceptorNumber ];
    electrodePositionsX = new double[electrodeNumber];
    electrodePositionsY = new double[electrodeNumber];

    occupation = new bool[acceptorNumber];
    for(int i=0;i<acceptorNumber;i++){
        occupation[i]=false;
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
        rates[i*(hoppingSiteNumber+1)]=0; //diagonal elements
        for(int j=acceptorNumber;j<i;j++){
            deltaEnergies[i*hoppingSiteNumber+j]=energies[j]-energies[i];
            deltaEnergies[j*hoppingSiteNumber+i]=energies[i]-energies[j];
            // std::cout<<"deltaEnergies el el "<< distances[i*hoppingSiteNumber+j]<<std::endl;

            if (deltaEnergies[i*hoppingSiteNumber+j] < 0){
                rates[i*hoppingSiteNumber+j]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA);
                rates[j*hoppingSiteNumber+i]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA-deltaEnergies[j*hoppingSiteNumber+i]);
            }
            else if (deltaEnergies[i*hoppingSiteNumber+j] > 0){
                rates[i*hoppingSiteNumber+j]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA-deltaEnergies[i*hoppingSiteNumber+j]);
                rates[j*hoppingSiteNumber+i]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA);
            }
            else{
                rates[i*hoppingSiteNumber+j]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA);
                rates[j*hoppingSiteNumber+i]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA);
            }
            constantRatesSumPart+=rates[i*hoppingSiteNumber+j];
            constantRatesSumPart+=rates[j*hoppingSiteNumber+i];
        }
    }

    readyForRun=true;

    DEBUG_FUNC_END
}


void System::updateRatesMPStoring(){ //same function as updateRatesSPStoring only mutex added
    DEBUG_FUNC_START

    mutex->lock_shared();
    if (knownRatesSum->count(hasedCurrentState)){ //state known
        ratesSum         = knownRatesSum->at(hasedCurrentState);
        partRatesSumList = konwnPartRatesSumList->at(hasedCurrentState);
        mutex->unlock_shared();

        ratesInMemory=true; //tell findSwap to use binary search

        // std::cout<<"found state: "<<hasedCurrentState<<std::endl;
    }
    else{ //state unknown
        mutex->unlock_shared();

        if (*storeKnownStates){ //still saving (memory limit not exceeded)
            //following part is nearly calcRates(), only ratesSum is not calculated (set as last component of partRatesSumList afterwards)
            {
                //acc acc hopp
                for(int i=0;i<acceptorNumber;i++){
                    if (occupation[i]){
                        for(int j=0;j<acceptorNumber;j++){
                            if (not occupation[j]){
                                // std::cout<<" ----- "<<i<<" "<<j<<"delta E "<<deltaEnergies[i*hoppingSiteNumber+j]<<std::endl;
                                deltaEnergies[i*hoppingSiteNumber+j]=energies[j]-energies[i]+pairEnergies[i*acceptorNumber+j];
                                // std::cout<<" ei "<<energies[i]<<" ej "<<energies[j]<<" epair "<<pairEnergies[i*acceptorNumber+j]<<std::endl;
                                if (deltaEnergies[i*hoppingSiteNumber+j] <= 0){
                                    rates[i*hoppingSiteNumber+j]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA);
                                }
                                else{
                                    rates[i*hoppingSiteNumber+j]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA-deltaEnergies[i*hoppingSiteNumber+j]);
                                }
                            }
                            else{
                                rates[i*hoppingSiteNumber+j]=0;
                                // deltaEnergies[i*hoppingSiteNumber+j]=0;
                            }
                        }
                    }
                    else{
                        for(int j=0;j<acceptorNumber;j++){ 
                            rates[i*hoppingSiteNumber+j]=0;
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
                                rates[i*hoppingSiteNumber+j]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA);
                            }
                            else{
                                rates[i*hoppingSiteNumber+j]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA-deltaEnergies[i*hoppingSiteNumber+j]);
                            }
                        }
                        else{
                            rates[i*hoppingSiteNumber+j]=0;
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
                                rates[i*hoppingSiteNumber+j]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA);
                            }
                            else{
                                rates[i*hoppingSiteNumber+j]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA-deltaEnergies[i*hoppingSiteNumber+j]);
                            }
                        }
                    }
                    else{
                        for(int j=acceptorNumber;j<hoppingSiteNumber;j++){ 
                            rates[i*hoppingSiteNumber+j]=0;
                            // deltaEnergies[i*hoppingSiteNumber+j]=0;
                        }
                    }
                }
            }

            partRatesSumList = std::make_shared<std::vector<double>>(hoppingSiteNumber*hoppingSiteNumber+1);
            (*partRatesSumList)[0]=0;
            for(int i=0; i < hoppingSiteNumber*hoppingSiteNumber; i++){
                (*partRatesSumList)[i+1]=(*partRatesSumList)[i]+rates[i];
            }
            ratesSum=(*partRatesSumList)[hoppingSiteNumber*hoppingSiteNumber];

            mutex->lock();
            (*knownRatesSum)[hasedCurrentState]=ratesSum;
            (*konwnPartRatesSumList)[hasedCurrentState]=partRatesSumList;
            mutex->unlock();
            ratesInMemory=true; //tell findSwap to use binary search


        }
        else{ //memory limit reached
            ratesInMemory=false; //tell findSwap to use comparison search
            updateRates();
        }
    }


    DEBUG_FUNC_END
}

void System::updateRatesSPStoring(){
    DEBUG_FUNC_START

    if (knownRatesSum->count(hasedCurrentState)){ //state known
        ratesSum         = knownRatesSum->at(hasedCurrentState);
        partRatesSumList = konwnPartRatesSumList->at(hasedCurrentState);
        ratesInMemory=true; //tell findSwap to use binary search

        // std::cout<<"found state: "<<hasedCurrentState<<std::endl;
    }
    else{ //state unknown
        if (*storeKnownStates){ //still saving (memory limit not exceeded)
            
            //following part is nearly calcRates(), only ratesSum is not calculated (set as last component of partRatesSumList afterwards)
            {
                //acc acc hopp
                for(int i=0;i<acceptorNumber;i++){
                    if (occupation[i]){
                        for(int j=0;j<acceptorNumber;j++){
                            if (not occupation[j]){
                                // std::cout<<" ----- "<<i<<" "<<j<<"delta E "<<deltaEnergies[i*hoppingSiteNumber+j]<<std::endl;
                                deltaEnergies[i*hoppingSiteNumber+j]=energies[j]-energies[i]+pairEnergies[i*acceptorNumber+j];
                                // std::cout<<" ei "<<energies[i]<<" ej "<<energies[j]<<" epair "<<pairEnergies[i*acceptorNumber+j]<<std::endl;
                                if (deltaEnergies[i*hoppingSiteNumber+j] <= 0){
                                    rates[i*hoppingSiteNumber+j]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA);
                                }
                                else{
                                    rates[i*hoppingSiteNumber+j]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA-deltaEnergies[i*hoppingSiteNumber+j]);
                                }
                            }
                            else{
                                rates[i*hoppingSiteNumber+j]=0;
                                // deltaEnergies[i*hoppingSiteNumber+j]=0;
                            }
                        }
                    }
                    else{
                        for(int j=0;j<acceptorNumber;j++){ 
                            rates[i*hoppingSiteNumber+j]=0;
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
                                rates[i*hoppingSiteNumber+j]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA);
                            }
                            else{
                                rates[i*hoppingSiteNumber+j]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA-deltaEnergies[i*hoppingSiteNumber+j]);
                            }
                        }
                        else{
                            rates[i*hoppingSiteNumber+j]=0;
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
                                rates[i*hoppingSiteNumber+j]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA);
                            }
                            else{
                                rates[i*hoppingSiteNumber+j]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA-deltaEnergies[i*hoppingSiteNumber+j]);
                            }
                        }
                    }
                    else{
                        for(int j=acceptorNumber;j<hoppingSiteNumber;j++){ 
                            rates[i*hoppingSiteNumber+j]=0;
                            // deltaEnergies[i*hoppingSiteNumber+j]=0;
                        }
                    }
                }
            }

            partRatesSumList = std::make_shared<std::vector<double>>(hoppingSiteNumber*hoppingSiteNumber+1);
            (*partRatesSumList)[0]=0;
            for(int i=0; i < hoppingSiteNumber*hoppingSiteNumber; i++){
                (*partRatesSumList)[i+1]=(*partRatesSumList)[i]+rates[i];
            }
            ratesSum=(*partRatesSumList)[hoppingSiteNumber*hoppingSiteNumber];

            (*knownRatesSum)[hasedCurrentState]=ratesSum;
            (*konwnPartRatesSumList)[hasedCurrentState]=partRatesSumList;
            ratesInMemory=true; //tell findSwap to use binary search

        }
        else{ //memory limit reached
            ratesInMemory=false; //tell findSwap to use comparison search
            updateRates();
        }
    }


    DEBUG_FUNC_END
}

void System::updateRates(){
    ratesSum = constantRatesSumPart;

    //acc acc hopp
    for(int i=0;i<acceptorNumber;i++){
        if (occupation[i]){
            for(int j=0;j<acceptorNumber;j++){
                if (not occupation[j]){
                    // std::cout<<" ----- "<<i<<" "<<j<<"delta E "<<deltaEnergies[i*hoppingSiteNumber+j]<<std::endl;
                    deltaEnergies[i*hoppingSiteNumber+j]=energies[j]-energies[i]+pairEnergies[i*acceptorNumber+j];
                    // std::cout<<" ei "<<energies[i]<<" ej "<<energies[j]<<" epair "<<pairEnergies[i*acceptorNumber+j]<<std::endl;
                    if (deltaEnergies[i*hoppingSiteNumber+j] <= 0){
                        rates[i*hoppingSiteNumber+j]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA);
                    }
                    else{
                        rates[i*hoppingSiteNumber+j]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA-deltaEnergies[i*hoppingSiteNumber+j]);
                    }
                    ratesSum+=rates[i*hoppingSiteNumber+j];
                }
                else{
                    rates[i*hoppingSiteNumber+j]=0;
                    // deltaEnergies[i*hoppingSiteNumber+j]=0;
                }
            }
        }
        else{
            for(int j=0;j<acceptorNumber;j++){ 
                rates[i*hoppingSiteNumber+j]=0;
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
                    rates[i*hoppingSiteNumber+j]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA);
                }
                else{
                    rates[i*hoppingSiteNumber+j]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA-deltaEnergies[i*hoppingSiteNumber+j]);
                }
                ratesSum+=rates[i*hoppingSiteNumber+j];
            }
            else{
                rates[i*hoppingSiteNumber+j]=0;
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
                    rates[i*hoppingSiteNumber+j]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA);
                }
                else{
                    rates[i*hoppingSiteNumber+j]=enhance::mediumFastExp(-2*distances[i*hoppingSiteNumber+j]/locLenA-deltaEnergies[i*hoppingSiteNumber+j]);
                }
                ratesSum+=rates[i*hoppingSiteNumber+j];
            }
        }
        else{
            for(int j=acceptorNumber;j<hoppingSiteNumber;j++){ 
                rates[i*hoppingSiteNumber+j]=0;
                // deltaEnergies[i*hoppingSiteNumber+j]=0;
            }
        }
    }
}


void System::updatePotential(){
    //this function needs to be splitted up in 3 parts for efficient parallel computing
    DEBUG_FUNC_START

    resetPotential();
    recalcPotential();
    setNewPotential();


    DEBUG_FUNC_END
}

void System::resetPotential(){
    DEBUG_FUNC_START
    
    //reset old potential
    for(int i=0;i<acceptorNumber;i++){
        energies[i]-=finEle->getPotential(acceptorPositionsX[i],acceptorPositionsY[i])*parameterStorage->parameters.at("e")/parameterStorage->parameters.at("kT");
    }
    for(int i=0;i<electrodeNumber;i++){
        energies[(i+acceptorNumber)]-=finEle->getPotential(electrodePositionsX[i],electrodePositionsY[i])*parameterStorage->parameters.at("e")/parameterStorage->parameters.at("kT");
    }

    DEBUG_FUNC_END
}

void System::recalcPotential(){
    DEBUG_FUNC_START

    //recalc potential
    for(int i=0;i < parameterStorage->electrodes.size();i++){
        finEle->updateElectrodeVoltage(i,parameterStorage->electrodes[i].voltage);
    }
    finEle->run();

    DEBUG_FUNC_END
}

void System::setNewPotential(){
    DEBUG_FUNC_START

    //set new potential
    for(int i=0;i<acceptorNumber;i++){
        energies[i]+=finEle->getPotential(acceptorPositionsX[i],acceptorPositionsY[i])*parameterStorage->parameters.at("e")/parameterStorage->parameters.at("kT");
    }
    for(int i=0;i<electrodeNumber;i++){
        energies[(i+acceptorNumber)]+=finEle->getPotential(electrodePositionsX[i],electrodePositionsY[i])*parameterStorage->parameters.at("e")/parameterStorage->parameters.at("kT");
    }


    // swapTrackFile.close(); // swapTracker
    // swapTrackFile.open(std::string("swapTrackFile")+std::to_string(fileNumber)+std::string(".txt"), ios::out); // swapTracker
    // fileNumber++; // swapTracker

    DEBUG_FUNC_END
}

void System::findSwap(){
    DEBUG_FUNC_START
    double rndNumber=enhance::random_double(0,ratesSum);
    double partRatesSum=0;
    for(int k=0; k<hoppingSiteNumber*hoppingSiteNumber;k++){
        partRatesSum+=rates[k];
        if(partRatesSum > rndNumber){
            lastSwapped1=k/hoppingSiteNumber;
            lastSwapped2=k%hoppingSiteNumber;

            std::cout<<"swapped1 "<<lastSwapped1<<" "<<lastSwapped2<<" k: "<<k<<std::endl;
            break;
        }
    }    
    // std::cout<<"rates: "<<std::endl;
    // for(int k=0; k<hoppingSiteNumber*hoppingSiteNumber;k++){
    //     std::cout<<rates[k]<<" ";
    // }
    // std::cout<<std::endl;
    // std::cout<<"ratesSum: "<<ratesSum<<" rndNumber: "<<rndNumber<<std::endl;

    DEBUG_FUNC_END
}

void System::findSwapBS(){
    DEBUG_FUNC_START
    double rndNumber=enhance::random_double(0,ratesSum);
    int l=0, r=hoppingSiteNumber*hoppingSiteNumber+1;
    int mid=-1;


    //binary search algorithm
    while (true){
        mid=(l+r)/2;
        // std::cout<<mid<<" ";

        if((*partRatesSumList)[mid]>rndNumber){
            if((*partRatesSumList)[mid-1]<rndNumber){
                lastSwapped1=(mid-1)/hoppingSiteNumber;
                lastSwapped2=(mid-1)%hoppingSiteNumber;
                // std::cout<<"swapped2 "<<lastSwapped1<<" "<<lastSwapped2<<" k: "<<mid-1<<std::endl;
                break;
            }
            else{
                r=mid;
            }
        }
        else if ((*partRatesSumList)[mid+1]>rndNumber){
            lastSwapped1=(mid)/hoppingSiteNumber;
            lastSwapped2=(mid)%hoppingSiteNumber;
            // std::cout<<"swapped2 "<<lastSwapped1<<" "<<lastSwapped2<<" k: "<<mid<<std::endl;
            break;
        }
        else{
            l=mid;
        }
    }

    // std::cout<<"rates: "<<std::endl;
    // for(int k=0; k<hoppingSiteNumber*hoppingSiteNumber;k++){
    //     std::cout<<rates[k]<<" ";
    // }
    // std::cout<<std::endl;
    // std::cout<<"partRatesSumList: "<<std::endl;
    // for(int k=0; k<hoppingSiteNumber*hoppingSiteNumber+1;k++){
    //     std::cout<<(*partRatesSumList)[k]<<" ";
    // }
    // std::cout<<std::endl;
    // std::cout<<"ratesSum: "<<ratesSum<<" rndNumber: "<<rndNumber<<std::endl;

    DEBUG_FUNC_END
}

void System::updateAfterSwap(){
    DEBUG_FUNC_START

    currentCounter[lastSwapped1]--;
    currentCounter[lastSwapped2]++;

    // swapTrackFile<<lastSwapped1<<";"<<lastSwapped2<<std::endl; // swapTracker

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




void System::run(int steps){
    DEBUG_FUNC_START


    if (storingMode) {
        if (parameterStorage->parameters.at("threads") > 1 and parameterStorage->parameters.at("shareMemory")){
            //multi processor mode (shared memory)
            for(int i=0; i<steps;i++){
                // check if memory limit is exceeded
                if (*storeKnownStates & (i%1000 ==0) & (((hoppingSiteNumber*hoppingSiteNumber+2)*8*knownRatesSum->size()) >= (parameterStorage->parameters.at("memoryLimit")*1e6))){
                    *storeKnownStates=false;
                    // std::cout<<"memory limit exceeded, stopping to store states"<<std::endl;
                }
                updateRatesMPStoring();
                if (ratesInMemory){
                    findSwapBS();
                }
                else{
                    findSwap();
                }
                updateAfterSwap();
                increaseTime();
            }
        }
        else {
            //single processor mode (also used when mp and no shared memory)
            for(int i=0; i<steps;i++){
                // check if memory limit is exceeded
                if (*storeKnownStates & (i%1000 ==0) & (((hoppingSiteNumber*hoppingSiteNumber+2)*8*knownRatesSum->size()) >= (parameterStorage->parameters.at("memoryLimit")*1e6))){
                    *storeKnownStates=false;
                    // std::cout<<"memory limit exceeded, stopping to store states"<<std::endl;
                }
                updateRatesSPStoring();
                if (ratesInMemory){
                    findSwapBS();
                }
                else{
                    findSwap();
                }
                updateAfterSwap();
                increaseTime();
            }
        }
    }
    else {
        for(int i=0; i<steps;i++){
            updateRates();
            findSwap();
            updateAfterSwap();
            increaseTime();
        }
    }
    
    // std::cout<<"done! curr: "<<" "<<currentCounter[int(parameterStorage->parameters.at("outputElectrode")+parameterStorage->parameters["acceptorNumber"])]/time <<std::endl;
    DEBUG_FUNC_END
}


void System::increaseTime(){
    DEBUG_FUNC_START

    time+=std::log(enhance::random_double(0,1))/(-1*ratesSum);

    DEBUG_FUNC_END
}
