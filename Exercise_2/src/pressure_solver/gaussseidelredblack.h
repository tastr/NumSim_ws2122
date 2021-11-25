#pragma once

#include <vector>
#include <array>
#include "pressuresolver.h"



class GaussSeidelRedBlack:
     public PressureSolver
{
protected: 
        double tol;
 
public:
      GaussSeidelRedBlack(Discretization& discretization_);
      ~GaussSeidelRedBlack();
      //FieldVariable calculateP();
      virtual void calculateP();
      void calculateBalckTiles();
      void calculateRedTiles();
     
    
    
    
};
