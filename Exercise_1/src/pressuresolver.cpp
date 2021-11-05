#include "pressuresolver.h"
#include <cassert>
#include "array2d.h"
#include "discretization.h"


PressureSolver::PressureSolver(Discretization discretization)
:discretization_(discretization)
{
    
}

PressureSolver::~PressureSolver()
{
}

   float PressureSolver::abs_(float number)
   {
 if (number>=0)
  {
    return number;
  }else
  {
    return -number;
    
  }
   }


    void PressureSolver::setPressureBoundaries()
    {int i_max = discretization_.pressure.size()[0], j_max = discretization_.pressure.size()[1];
     
     for (int j = 0; j < j_max; j++)
     {
         discretization_.pressure(0,j)    = discretization_.pressure(1,j);
         discretization_.pressure(i_max-1,j)= discretization_.pressure(i_max-2,j);
     }
        
     for (int i = 0; i < i_max; i++)
     {
        discretization_.pressure(i,0)= discretization_.pressure(i,1);
        discretization_.pressure(i,j_max-1)= discretization_.pressure(i,j_max-2);
     }
     
    }