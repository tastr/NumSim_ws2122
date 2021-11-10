#include "centraldifferences.h"

CentralDifferences::CentralDifferences(Settings settings):Discretization(settings)
  { 
  }

 CentralDifferences::~CentralDifferences()
  {
  }


double CentralDifferences::computeDu2Dx(int i, int j) const    
 {double uiphj=(u(i,j)  + u(i+1,j)) ;
  double uimhj=u(i-1,j) + u(i,j);
   return ( uiphj*uiphj - uimhj*uimhj)/(4*delta_x); 
 } 

 double CentralDifferences::computeDv2Dy(int i, int j) const
 {double vijph = (v(i,j) + v(i,j+1));
  double vijmh = (v(i,j-1) + v(i,j));
   return ( vijph*vijph-vijmh*vijmh )/(4*delta_y);
}

 double CentralDifferences::computeDuvDx(int i, int j) const
 {double viphj = (v(i,j) + v(i,j+1)) * (u(i,j)  + u(i+1,j)) ;
  double vimhj  = (v(i-1,j) + v(i-1,j+1)) * (u(i-1,j)+ u(i,j));
    return (viphj - vimhj )/(delta_x*4);
 }

 double CentralDifferences::computeDuvDy(int i, int j) const
 {double viphj   = (v(i,j)+v(i+1,j)) * (u(i,j)+u(i,j+1));
  double viphjm1 = (v(i,j-1)+v(i+1,j-1)) * (u(i,j-1)+u(i,j));
   return (viphj-viphjm1)/(4*delta_y);
 }


// the index was fixed here as well
void CentralDifferences::calculation()
{
    for (int j = 1; j < settings_.nCells[1]+1; j++)
    {
        for (int i = 1; i < settings_.nCells[0]+1; i++)
        {   
            F(i,j) = u(i,j) + deltat * ((computeDuDx2(i,j) + computeDuDy2(i,j)) / settings_.re - computeDuvDy(i,j) - computeDu2Dx(i,j) + settings_.g[0]);
            G(i,j) = v(i,j) + deltat * ((computeDvDx2(i,j) + computeDvDy2(i,j)) / settings_.re - computeDuvDx(i,j) - computeDv2Dy(i,j) + settings_.g[1]);
        }
    }   
}