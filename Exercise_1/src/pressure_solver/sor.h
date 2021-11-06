#pragma once

#include <vector>
#include <array>
#include "pressuresolver.h"



class SOR:
     public PressureSolver
{
protected: 
         double omega; //Omega festgelegt durch parameter
 
public:
      SOR(Discretization discretization_);
      ~SOR();
      //FieldVariable Iterationsverfahren();
      void Iterationsverfahren();
      float tol;
     
    
    
    
};
