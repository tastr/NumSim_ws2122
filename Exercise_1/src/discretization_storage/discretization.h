#pragma once

#include "staggeredgrid.h"
#include "settings.h"
#include "fieldvariable.h"

class Discretization : 
    public StaggeredGrid
{
protected:
   double deltat;
   FieldVariable F, G, rhs_;
public:
    Discretization(Settings settings);
    virtual ~Discretization();
   void calculation_altfinitedifferenzen();

   virtual void calculation();
   void updateVelocity();
   void updateBoundaryFG(); //for driven cavity not neccessary, but for other cases maybe.


   void updateDeltaT();
   double  min2(double value1, double value2) const; 
   double  min3(double value1, double value2, double value3) const;
   
   double computeDuDx2(int i, int j) const;
   double computeDvDy2(int i, int j) const;
   double computeDuDy2(int i, int j) const;
   double computeDvDx2(int i, int j) const;

   double computeDuDx(int i, int j) const;
   double computeDvDy(int i, int j) const;

   double computeDpDx(int i, int j) const;
   double computeDpDy(int i, int j) const;

   // die koenten wir eigentlich entfernen
   virtual double computeDu2Dx(int i, int j) const ; 
   virtual double computeDv2Dy(int i, int j) const ;
   virtual double computeDuvDx(int i, int j) const ;
   virtual double computeDuvDy(int i, int j) const ;

   // get functions:
   
    double f(int i, int j) const; 
    double g(int i, int j) const; 
    double rhs(int i, int j) const; 
    double getOmega() const;
    double getDeltaT() const;

    // set functions
    void setRHS(int i, int j, double value);
};


