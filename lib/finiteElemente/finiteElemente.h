// TODO: update voltages and not replace whole system


#ifndef FINITE_ELEMENTE_H
#define FINITE_ELEMENTE_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <numeric>

#include "mfem.hpp"

#define PI 3.14159265

using namespace mfem;

class FiniteElementeBase
{

protected:
    bool saveSolution;
    int runNumber=0;
    int const dim = 2, sdim = 2, order =1;
    std::vector<std::vector<int>> electrodeVertexIndices;
    int numberOfElectrodes=0;

    Mesh *mesh;

    //mfem stuff
    FiniteElementCollection *fec;
    FiniteElementSpace *fespace;
    BilinearForm * a; 
    LinearForm * b;
    OperatorPtr A;
    Vector B, X;
    GSSmoother M;
    Array<int> ess_tdof_list;


public:
    FiniteElementeBase(bool saveSolution);
    ~FiniteElementeBase();
    GridFunction * solutionVector; // changed to pointer, named x in example

    void initRun(bool initDevice = false);
    void run();
    void updateElectrodeVoltage(int const & electrodeIndex, double const & voltage);

    virtual double getPotential(double const & x, double const & y){};
    virtual void setElectrode(double const & voltage, double begin, double end){}; //begin/end given in radiant
    virtual void setElectrode(double const & voltage, double begin, double end, int edge ){}; //begin/end given in length units

};



class FiniteElementeRect : public FiniteElementeBase
{

private:
    double const len = 0, width = 0;
    int numberVerticesX,numberVerticesY;
    int** vertexIndexMap;
    void initMesh  (int const & maxNumberOfElements);

public:
    ~FiniteElementeRect();
    FiniteElementeRect(double const & len,   double const & width, int const & maxNumberOfElments, bool saveSolution=false); 

    void setElectrode(double const & voltage, double begin, double end, int edge ); //begin/end given in length units
    double getPotential(double const & x, double const & y); //nearest neighbour interpolation
};



class FiniteElementeCircle : public FiniteElementeBase
{

private:
    double const radius = 0;
    int layers; // only used in circ mode
    void initMesh(int const & maxNumberOfElements);

public:
    FiniteElementeCircle(double const & radius, int const & maxNumberOfElments, bool saveSolution=false);

    void setElectrode(double const & voltage, double begin, double end); //begin/end given in radiant
    double getPotential(double const & x, double const & y);
};


#endif // FINITE_ELEMENTE_H
