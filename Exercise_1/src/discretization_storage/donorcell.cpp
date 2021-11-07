#include "donorcell.h"
#include "discretization.h"

   DonorCell::DonorCell(Settings settings):Discretization(settings) 
   {
   }

   DonorCell::~DonorCell() 
   {
   }

double DonorCell::computeDu2Dx(int i, int j) const    
 {double uiphj      = u(i,j)       + u(i+1,j);
  double uimhj      = u(i-1,j)     + u(i,j);
  double donorTerm1 = abs(u(i,j)   + u(i+1,j)) * (u(i,j)-u(i+1,j));
  double donorTerm2 = abs(u(i-1,j) + u(i,j))   * (u(i-1,j)-u(i,j));
  std::cout << "Donorcell" <<std::endl;
   return ( uiphj*uiphj + uiphj*uiphj + settings_.alpha * ( donorTerm1- donorTerm2) )/(4*delta_x); 
 } 

 double DonorCell::computeDv2Dy(int i, int j) const
 {double uijph = (v(i,j)   + v(i,j+1));
  double uijmh = (v(i,j-1) + v(i,j));
  double donorTerm1 = abs(v(i,j) + v(i,j+1)) * (v(i,j)-v(i,j+1));
  double donorTerm2 = abs(v(i,j-1) + v(i,j)) * (v(i,j-1)-v(i,j));
  std::cout << "Donorcell" <<std::endl;
   return ( uijph-uijmh + settings_.alpha * (donorTerm1 - donorTerm2) )/(4*delta_x);
 }

 double DonorCell::computeDuvDx(int i, int j) const
 {double donorTerm1 = (v(i+1,j)-v(i,j))*abs(u(i,j+1)+u(i,j));
  double donorTerm2 = (v(i,j)-v(i-1,j))*abs(u(i-1,j+1) + u(i-1,j)); 
  double viphjm1  = v(i+1,j-1) + v(i,j-1) ;
  double viphj    = v(i+1,j)   + v(i,j)   ;
  double uijph = (v(i,j) + v(i,j+1));
  double uijmh = (v(i,j-1) + v(i,j));
   return ((viphj * uijph - viphjm1 * uijmh) + settings_.alpha *(donorTerm1-donorTerm2))/(delta_y*4);
 }

 double DonorCell::computeDuvDy(int i, int j) const
 {double viphj = (v(i,j)+v(i+1,j));
  double viphjm1 = (v(i,j-1)+v(i+1,j-1));
  double donorTerm1 = abs(v(i,j)+v(i+1,j))     * (u(i,j)-u(i,j+1));
  double donorTerm2 = abs(v(i,j-1)+v(i+1,j-1)) * (u(i,j-1)-u(i,j)); 
   return (viphj-viphjm1 + settings_.alpha * (donorTerm1-donorTerm2))/(4*delta_y);
 }



void DonorCell::calculation()
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


