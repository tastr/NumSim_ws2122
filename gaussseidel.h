#pragma once

#include <vector>
#include <array>
#include "pressuresolver.h"



class GaussSeidel:
     public PressureSolver
{
protected: 
	 Discretization discretization;
        double tol;
 
public:
      GaussSeidel(Discretization discretization_);
      ~GausSeidel();
      //FieldVariable Iterationsverfahren();
      void Iterationsverfahren();
     
     
    
    
    
};
