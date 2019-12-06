#include <iostream>
#include <fstream>

#include <vector>
#include <limits.h>
#include "finiteElemente.h"


int main(int argc, char *argv[]){

    FiniteElemente finEle(1,1,1e2);


    // finEle.setElectrode(0,0.5,0,10);
    // finEle.setElectrode(0,0.5,2,10);
    // // finEle.setElectrode(0.4,0.5,1,-10);
    // finEle.setElectrode(0.2,0.3,2,-10);
    // // finEle.setElectrode(0.7,0.8,3,10);


    finEle.setElectrode(0.2,0.3,0, 1);
    finEle.setElectrode(0.7,0.8,0,-1);
    finEle.setElectrode(0.2,0.3,1, 1);
    finEle.setElectrode(0.7,0.8,1,-1);
    finEle.setElectrode(0.2,0.3,2, 1);
    finEle.setElectrode(0.7,0.8,2,-1);
    finEle.setElectrode(0.2,0.3,3, 1);
    finEle.setElectrode(0.7,0.8,3,-1);

    finEle.initRun();
    finEle.run();

    finEle.updateElectrodeVoltage(1,-10);
    finEle.updateElectrodeVoltage(3,-10);

    finEle.run();




    // std::ofstream f("out.dat",std::ofstream::trunc);
    // int res=1000;
    // for(int i=0;i<res;i++){
    //     for(int j=0;j<res;j++){
            
    //     }
    // }


    // std::cout<<"pot "<<finEle.getPotential(0.7,0.2)<<std::endl;



   return 0;
}
