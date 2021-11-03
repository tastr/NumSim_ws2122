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
};


