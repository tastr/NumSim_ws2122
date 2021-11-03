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
};


