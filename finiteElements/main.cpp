#include <iostream>
#include <fstream>

#include <vector>

#include "finiteElemente.h"


int main(int argc, char *argv[]){

    FiniteElemente finEle(1,1,10);
    finEle.setElectrode(0.1,0.9,0,10);


    finEle.printMesh("test.mesh");

    // finEle.run();



   return 0;
}
