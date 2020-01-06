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
    bool saveSolution;
    int runNumber=0;
    int numberVerticesX,numberVerticesY;
    int const dim = 2, sdim = 2, order =1;
    int** vertexIndexMap;
    std::vector<std::vector<int>> electrodeVertexIndices;
    int numberOfElectrodes=0;
    void initMesh(int const & maxNumberOfElements);

    Mesh *mesh;

    //mfem stuff
    FiniteElementCollection *fec;
    FiniteElementSpace *fespace;
    GridFunction *solutionVector; // changed to pointer, named x in example
    BilinearForm * a; 
    LinearForm * b;
    OperatorPtr A;
    Vector B, X;
    GSSmoother M;
   Array<int> ess_tdof_list;


public:
    FiniteElemente(double const & len,double const & width, int const & maxNumberOfElments, bool saveSolution=false);
    ~FiniteElemente();

    void initRun();
    void run();
    void setElectrode(double const & begin, double const & end, int const & edge, double const & voltage);
    void updateElectrodeVoltage(int const & electrodeIndex, double const & voltage);
    double getPotential(double const & x, double const & y); //nearest neighbour interpolation
   
};

#endif // FINITE_ELEMENTE_H
