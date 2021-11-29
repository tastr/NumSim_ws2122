#pragma once

#include "staggeredgrid.h"
#include "../settings.h"
#include "fieldvariable.h"
#include "../partitioning/partitioning.h"
#include <mpi.h>


class Discretization : 
    public StaggeredGrid
{
protected:
   double deltat;
   FieldVariable F, G, rhs_;
//    Partitioning partitioning_;
public:
   Discretization(Settings settings, Partitioning partitioning);
   virtual ~Discretization();
   void calculation_altfinitedifferenzen();

   virtual void calculation() = 0;
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
   virtual double computeDu2Dx(int i, int j) const = 0; 
   virtual double computeDv2Dy(int i, int j) const = 0; 
   virtual double computeDuvDx(int i, int j) const = 0; 
   virtual double computeDuvDy(int i, int j) const = 0; 

   // get functions:
   
    double f(int i, int j) const; 
    double g(int i, int j) const; 
    double rhs(int i, int j) const; 
    double getOmega() const;
    double getDeltaT() const;

    // set functions
    void setRHS(int i, int j, double value);

    void setPressureBCParalell();
    void setBorderVelocityParalell(std::array<double,2> top,std::array<double,2> left,std::array<double,2> right,std::array<double,2> bottom);
    void updateBoundaryFGParalell(); //for driven cavity not neccessary, but for other cases maybe.

};


