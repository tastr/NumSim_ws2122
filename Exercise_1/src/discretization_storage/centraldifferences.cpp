#include "centraldifferences.h"

CentralDifferences::CentralDifferences(Settings settings):Discretization(settings)
  {
    
  }

 CentralDifferences::~CentralDifferences()
  {
  }


double CentralDifferences::computeDu2Dx(int i, int j) const    
 {double uiphj=(u(i,j) + u(i+1,j)) ;
  double uimhj=u(i-1,j) + u(i,j);
   return ( uiphj*uiphj + uiphj*uiphj)/(4*delta_x); 
 } 

 double CentralDifferences::computeDv2Dy(int i, int j) const
 {double uijph = (v(i,j) + v(i,j+1));
  double uijmh = (v(i,j-1) + v(i,j));
  std::cout << "Zentrale Differenzen" <<std::endl;
   return ( uijph-uijmh )/(4*delta_x);
}

 double CentralDifferences::computeDuvDx(int i, int j) const
 {double viphjm1  = v(i+1,j-1) + v(i,j-1) ;
  double viphj    = v(i+1,j)   + v(i,j)   ;
  double uijph = (v(i,j) + v(i,j+1));
  double uijmh = (v(i,j-1) + v(i,j));
  std::cout << "Zentrale differenzen" <<std::endl;
   return ((viphj * uijph - viphjm1 * uijmh) )/(delta_y*4);
 }

 double CentralDifferences::computeDuvDy(int i, int j) const
 {double viphj = (v(i,j)+v(i+1,j));
  double viphjm1 = (v(i,j-1)+v(i+1,j-1));
   return (viphj-viphjm1)/(4*delta_y);
 }



void CentralDifferences::calculation()
{
    for (int j = 1; j < settings_.nCells[1]; j++)
    {
        for (int i = 1; i < settings_.nCells[0]; i++)
        {   
            F(i,j) = u(i,j) + deltat * ((computeDuDx2(i,j) + computeDuDy2(i,j)) / settings_.re - computeDuvDy(i,j) - computeDu2Dx(i,j) + settings_.g[0]);
            G(i,j) = v(i,j) + deltat * ((computeDvDx2(i,j) + computeDvDy2(i,j)) / settings_.re - computeDuvDx(i,j) - computeDv2Dy(i,j) + settings_.g[1]);
        }
    }   
}