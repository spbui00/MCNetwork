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


/*!
  Class to solve laplace equation.
  Parent class istself cant run, only childs.
  Analogously implemented to https://github.com/mfem/mfem/blob/master/examples/ex1.cpp.
*/
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
    /*!
        \param saveSolution if true: save mesh and solution in each step. can be visualized using "glivs -m finEle.mesh -g laplace_solution0.gf"
    */
    FiniteElementeBase(bool saveSolution);
    ~FiniteElementeBase();
    GridFunction * solutionVector; // changed to pointer, named x in example

    void initRun(bool initDevice = false);
    void run();
    void updateElectrodeVoltage(int const & electrodeIndex, double const & voltage);
    /*!
        Get Potential using nearest neighbour interpolation
        \param x in length units
        \param y in length units
    */
    virtual double getPotential(double const & x, double const & y){};
    /*!
        Set electrode position and voltage (boundary condition of laplace equation), polar version.
        \param voltage in volts
        \param begin in radiant
        \param end in radiant
    */
    virtual void setElectrode(double const & voltage, double begin, double end){};
    /*!
        Set electrode position and voltage (boundary condition of laplace equation), cartesian version.
        \param voltage in volts
        \param begin in length units
        \param end in length units
        \param edge see ElectrodeParameters::edge
    */
   virtual void setElectrode(double const & voltage, double begin, double end, int edge ){};

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
    double deltaR;
    int layers; // only used in circ mode
    void initMesh(int const & maxNumberOfElements);

public:
    FiniteElementeCircle(double const & radius, int const & maxNumberOfElments, bool saveSolution=false);

    void setElectrode(double const & voltage, double begin, double end); //begin/end given in radiant
    double getPotential(double const & x, double const & y);
};


#endif // FINITE_ELEMENTE_H
