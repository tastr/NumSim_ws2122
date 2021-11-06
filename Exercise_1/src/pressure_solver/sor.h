#pragma once

#include <vector>
#include <array>
#include "pressuresolver.h"



class SOR:
     public PressureSolver
{
protected: 
	 Discretization discretization;
         double omega;
 
public:
      SOR(Discretization discretization_);
      ~SOR();
      //FieldVariable Iterationsverfahren();
      void Iterationsverfahren();
      float tol;
     
    
    
    
};
