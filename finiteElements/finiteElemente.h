#ifndef FINITE_ELEMENTE_H
#define FINITE_ELEMENTE_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include "mfem.hpp"

using namespace std;
using namespace mfem;

class FiniteElemente
{

private:
    double len, width;
    int dim = 2, sdim = 2, order =1, attrib = 1; 
    int numberVerticesX,numberVerticesY;

    Mesh *mesh;

public:
    FiniteElemente(double len_,double width_, int maxNumberOfElments);
    ~FiniteElemente();
    void printMesh(std::string);
    void run();
    void setElectrode(double a, double b, int edge,double voltage);

   
};

#endif // FINITE_ELEMENTE_H
