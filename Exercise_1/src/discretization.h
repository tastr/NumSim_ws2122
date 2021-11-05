#pragma once

#include "staggeredgrid.h"
#include "settings.h"
#include "fieldvariable.h"

class Discretization : 
    public StaggeredGrid
{
private:
   // StaggeredGrid theGrid;
   Settings settings_;

public:
    Discretization(Settings settings);
    ~Discretization();
   void calculation();
   float min2(float value1, float value2) const; 
   float min3(float value1, float value2, float value3) const;
<<<<<<< HEAD
   
   float deltat;
   float delta_x;
   float delta_y;
=======

   // get functions:
   // oder doch in staggered grid definieren?
    float f(int i, int j) const; //TODO
    float g(int i, int j) const; //TODO
    float rhs(int i, int j) const; //TODO
>>>>>>> b58dd6ee07f3c0424f80da1333cd28f86fbe51bf
};


