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
    virtual ~Discretization();
   void calculation_altfinitedifferenzen();

   virtual void calculation();
   
   void updateDeltaT();
   double  min2(double value1, double value2) const; 
   double  min3(double value1, double value2, double value3) const;
   
   double computeDuDx2(int i, int j) const;
   double computeDvDy2(int i, int j) const;
   double computeDuDy2(int i, int j) const;
   double computeDvDx2(int i, int j) const;

   double computeDuDx(int i, int j) const;
   double computeDvDy(int i, int j) const;

   virtual double computeDu2Dx(int i, int j) const ; // keine Ahnung ob das = 0 hier hin gehoert
   virtual double computeDv2Dy(int i, int j) const ;
   virtual double computeDuvDx(int i, int j) const ;
   virtual double computeDuvDy(int i, int j) const ;

   // get functions:
   // oder doch in staggered grid definieren?
    double f(int i, int j) const; 
    double g(int i, int j) const; 
    double rhs(int i, int j) const; //const statment removed becoause of error
    double getOmega() const; //TODO
    double getDeltaT() const;

    // set functions
    void setRHS(int i, int j, double value);
};


