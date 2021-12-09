#pragma once

#include <vector>
#include <array>
#include "pressuresolver.h"



class GaussSeidel:
     public PressureSolver
{
protected: 
        double tol;
 
public:
      GaussSeidel(Discretization& discretization_);
      ~GaussSeidel();
      //FieldVariable calculateP();
      virtual double calculateP();
     
     
    
    
    
};
