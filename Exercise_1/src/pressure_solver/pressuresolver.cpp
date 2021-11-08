#include "pressuresolver.h"
#include <cassert>

PressureSolver::PressureSolver(Discretization& discretization)
:discretization_(discretization)
{
    
}

PressureSolver::~PressureSolver()
{
}

double PressureSolver::abs_(double number)
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
    {int i_max = discretization_.getSize()[0], j_max = discretization_.getSize()[1];
     
     for (int j = 0; j < j_max; j++)
     {
         discretization_.setP(0,j,discretization_.p(1,j));
         discretization_.setP(i_max-1,j,discretization_.p(i_max-2,j));
     }
        
     for (int i = 0; i < i_max; i++)
     {
        discretization_.setP(i,0,discretization_.p(i,1));
        discretization_.setP(i,j_max-1,discretization_.p(i,j_max-2));
     }
    }


void PressureSolver::calculateRHS()
{
  double current_rhs=0;
  for (int j = 1; j < discretization_.getSize()[1]-1; j++)
  {
    for (int i = 1; i < discretization_.getSize()[0]-1; i++)
    {
      current_rhs=((discretization_.f(i,j)-discretization_.f(i-1,j))/discretization_.dx()+(discretization_.g(i,j)-discretization_.g(i,j-1))/discretization_.dy())/discretization_.getDeltaT();
      discretization_.setRHS(i,j, current_rhs);
    }
  }
}

void PressureSolver::calculateP()
{
  //is virtual class only should never be called
  assert(false);
}

double PressureSolver::residuum() 
{ double res =0;
  for (int j = 1; j < discretization_.getSize()[1]-1; j++)
  {
      for (int i = 1; i < discretization_.getSize()[0]-1; i++)
      {
       res = res + abs_( (discretization_.p(i-1,j)  - 2 * discretization_.p(i,j) + discretization_.p(i+1,j))/(discretization_.dx() * discretization_.dx()) + (discretization_.p(i,j-1) - 2 * discretization_.p(i,j)  + discretization_.p(i,j+1)) / (discretization_.dy() * discretization_.dy()) - discretization_.rhs(i,j));
      }     
  }
  return res;
}


