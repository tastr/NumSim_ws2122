#pragma once

#include "staggeredgrid.h"
#include "settings.h"
#include "fieldvariable.h"

class Discretization : 
    public StaggeredGrid
{
protected:
   // StaggeredGrid theGrid;
//    Settings settings_;
   double deltat;
   FieldVariable F, G, rhs_;
public:
    Discretization(Settings settings);
    ~Discretization();
   void calculation();
   void updateDeltaT();
   double  min2(double value1, double value2) const; 
   double  min3(double value1, double value2, double value3) const;
   

   // get functions:
   // oder doch in staggered grid definieren?
    double f(int i, int j) const; 
    double g(int i, int j) const; 
    double rhs(int i, int j) const; //const statment removed becoause of error
    double getOmega() const; //TODO
    double getDeltaT() const;
};


