#pragma once

#include <vector>
#include <array>
#include "pressuresolver.h"



class SOR:
     public PressureSolver
{
protected: 
	 Discretization discretization;
 
 
public:
      SOR(Discretization discretization_);
      ~SOR();
      //FieldVariable Iterationsverfahren();
      void Iterationsverfahren();
      float tol;
     
    
    
    
};
