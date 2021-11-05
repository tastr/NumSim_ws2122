#pragma once

#include <vector>
#include <array>
#include "discretization.h"


class PressureSolver
{
protected: 
     
 
 
public:
      PressureSolver(Discretization discretization);
      ~PressureSolver();
       float abs_(float number);
       void setPressureBoundaries();
       Discretization discretization_;
       
};
