#include <iostream>
#include <fstream>

#include <vector>
#include <limits.h>
#include "finiteElemente.h"
#include <math.h>
#include <random>


int main(int argc, char *argv[]){



    FiniteElementeCircle finEle(300, 5e3, true);


    double electrodeWidth = 60/300.0;
    int i;


    i = 0;
    finEle.setElectrode(1  , i*2*PI/8-electrodeWidth/2, i*2*PI/8+electrodeWidth/2);
    i = 1;
    finEle.setElectrode(-1 , i*2*PI/8-electrodeWidth/2, i*2*PI/8+electrodeWidth/2);
    i = 2;
    finEle.setElectrode(0.3, i*2*PI/8-electrodeWidth/2, i*2*PI/8+electrodeWidth/2);
    i = 3;
    finEle.setElectrode(-0.2, i*2*PI/8-electrodeWidth/2, i*2*PI/8+electrodeWidth/2);
    i = 4;
    finEle.setElectrode(0.8, i*2*PI/8-electrodeWidth/2, i*2*PI/8+electrodeWidth/2);
    i = 5;
    finEle.setElectrode(0.7, i*2*PI/8-electrodeWidth/2, i*2*PI/8+electrodeWidth/2);
    i = 6;
    finEle.setElectrode(-1, i*2*PI/8-electrodeWidth/2, i*2*PI/8+electrodeWidth/2);
    i = 7;
    finEle.setElectrode(0, i*2*PI/8-electrodeWidth/2, i*2*PI/8+electrodeWidth/2);

    finEle.initRun();
    finEle.run();



    // FiniteElementeRect finEle(1,1,1e2, true);



    // finEle.setElectrode( 1,0.2,0.3,0);
    // finEle.setElectrode(-1,0.7,0.8,0);
    // // finEle.setElectrode( 1,0.2,0.3,1);
    // // finEle.setElectrode(-1,0.7,0.8,1);
    // // finEle.setElectrode( 1,0.2,0.3,2);
    // // finEle.setElectrode(-1,0.7,0.8,2);
    // // finEle.setElectrode( 1,0.2,0.3,3);
    // // finEle.setElectrode(-1,0.7,0.8,3);

    // finEle.initRun();
    // finEle.run();

    // finEle.updateElectrodeVoltage(1,-10);
    // finEle.updateElectrodeVoltage(3,-10);

    // finEle.run();




    // std::ofstream f("out.dat",std::ofstream::trunc);
    // int res=1000;
    // for(int i=0;i<res;i++){
    //     for(int j=0;j<res;j++){
            
    //     }
    // }


    // std::cout<<"pot "<<finEle.getPotential(0.7,0.2)<<std::endl;



   return 0;
}
