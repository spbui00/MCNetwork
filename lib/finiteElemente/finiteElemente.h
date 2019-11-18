// TODO: update voltages and not replace whole system


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
    int numberVerticesX,numberVerticesY;
    int dim = 2, sdim = 2, order =1;
    int** vertexIndexMap;
    void initMesh(int maxNumberOfElements);

    Mesh *mesh;

    FiniteElementCollection *fec;
    FiniteElementSpace *fespace;
    GridFunction *solutionVector; // changed to pointer

public:
    FiniteElemente(double len_,double width_, int maxNumberOfElments);
    ~FiniteElemente();

    void run();
    void setElectrode(double a, double b, int edge,double voltage);
    double getPotential(double x, double y); //nearest neighbour interpolation
   
};

#endif // FINITE_ELEMENTE_H
