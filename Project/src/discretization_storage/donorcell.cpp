#include "donorcell.h"
//#include "discretization.h"

   DonorCell::DonorCell(Settings settings):Discretization(settings) 
   {
   }

   DonorCell::~DonorCell() 
   {
   }

//functions for numerical derivation using the donorcell scheme 
double DonorCell::computeDu2Dx(int i, int j) const    
 {double uiphj=(u(i,j)  + u(i+1,j)) ;
  double uimhj=u(i-1,j) + u(i,j);
  double donorTerm1 = abs(u(i,j)   + u(i+1,j)) * (u(i,j)-u(i+1,j));
  double donorTerm2 = abs(u(i-1,j) + u(i,j))   * (u(i-1,j)-u(i,j));
  return  (uiphj*uiphj - uimhj*uimhj + settings_.alpha * ( donorTerm1- donorTerm2) )/(4*delta_x); 
  

 
 
  //double uiphj=(u(i,j)  + u(i+1,j)) ;
  //double uimhj=u(i-1,j) + u(i,j);
  //double donorTerm1 = abs(u(i,j)   + u(i+1,j)) * (u(i,j)-u(i+1,j));
  //double donorTerm2 = abs(u(i-1,j) + u(i,j))   * (u(i-1,j)-u(i,j));
   
    //return ( uiphj*uiphj - uimhj*uimhj + settings_.alpha * ( donorTerm1- donorTerm2) )/(4*delta_x); 
 } 

 double DonorCell::computeDv2Dy(int i, int j) const
 {double vijph = (v(i,j) + v(i,j+1));
  double vijmh = (v(i,j-1) + v(i,j));
//return ( vijph*vijph - vijmh*vijmh )/(4*delta_y);
  double donorTerm1 = abs(v(i,j) + v(i,j+1)) * (v(i,j)-v(i,j+1));
  double donorTerm2 = abs(v(i,j-1) + v(i,j)) * (v(i,j-1)-v(i,j));
  return ( vijph*vijph-vijmh*vijmh + settings_.alpha * (donorTerm1 - donorTerm2) )/(4*delta_y);
  
     
  //double vijph = (v(i,j)   + v(i,j+1));
  //double vijmh = (v(i,j-1) + v(i,j));
  //double donorTerm1 = abs(v(i,j) + v(i,j+1)) * (v(i,j)-v(i,j+1));
  //double donorTerm2 = abs(v(i,j-1) + v(i,j)) * (v(i,j-1)-v(i,j));
   
   //return ( vijph*vijph-vijmh*vijmh + settings_.alpha * (donorTerm1 - donorTerm2) )/(4*delta_y);
 }

 double DonorCell::computeDuvDx(int i, int j) const
 {  double v_iphj     = (v(i+1,j)+v(i,j));
    double u_ijph     = (u(i,j+1)+u(i,j));
    double v_imhj     = (v(i,j)+v(i-1,j));
    double u_imhjph   = (u(i-1,j+1)+u(i-1,j));
    double result_neu = (v_iphj*u_ijph) - (v_imhj*u_imhjph);
    double donorTerm1 = (v(i+1,j)-v(i,j))*abs(u_ijph);
    double donorTerm2 = (v(i,j)-v(i-1,j))*abs(u_imhjph); 
    
    return (result_neu  - settings_.alpha *(donorTerm1-donorTerm2))/(delta_x*4); //nicht sicher, warum mit minus besser....
 
  
  //double donorTerm1 = (v(i,j)-v(i,j+1))*abs(u(i,j)+u(i+1,j));
  //double donorTerm2 = (v(i-1,j)-v(i-1,j+1))*abs(u(i-1,j) + u(i,j)); 
  //double viphj  = (v(i,j) + v(i,j+1)) * (u(i,j)  + u(i+1,j)) ;
  //double vimhj  = (v(i-1,j) + v(i-1,j+1)) * (u(i-1,j)+ u(i,j));
   
   //return ((viphj - vimhj ) + settings_.alpha *(donorTerm1-donorTerm2))/(delta_x*4);
 }

 double DonorCell::computeDuvDy(int i, int j) const
 {double viphj   = (v(i,j)+v(i+1,j)) * (u(i,j)+u(i,j+1));
  double viphjm1 = (v(i,j-1)+v(i+1,j-1)) * (u(i,j-1)+u(i,j));
  double donorTerm1  = abs(v(i,j)+v(i+1,j))     * (u(i,j)-u(i,j+1));
  double donorTerm2  = abs(v(i,j-1)+v(i+1,j-1)) * (u(i,j-1)-u(i,j)); 
  
  return (viphj-viphjm1 + settings_.alpha * (donorTerm1-donorTerm2))/(4*delta_y);
   
    
    //return (viphj-viphjm1)/(4*delta_y);
   
  //double viphj       = (v(i,j)+v(i+1,j))        * (u(i,j)+u(i,j+1));
  //double viphjm1     = (v(i,j-1)+v(i+1,j-1))    * (u(i,j-1)+u(i,j));
  //double donorTerm1  = abs(v(i,j)+v(i+1,j))     * (u(i,j)-u(i,j+1));
  //double donorTerm2  = abs(v(i,j-1)+v(i+1,j-1)) * (u(i,j-1)-u(i,j)); 
   
   //return (viphj-viphjm1 + settings_.alpha * (donorTerm1-donorTerm2))/(4*delta_y);
 }


// calculates the values of auxilliarz variables F and G using the specific numerical
// derivative functions of the class


// #TODO
void DonorCell::calculation()
{
    for (int j = uJBegin()+1; j < uJEnd()-1; j++)
    {
        for (int i = uIBegin()+1; i < uIEnd()-1; i++)
        {
          if (type(i,j)==0)
          {  
            F(i,j) = u(i,j) + deltat * ((computeDuDx2(i,j) + computeDuDy2(i,j)) / settings_.re - computeDuvDy(i,j) - computeDu2Dx(i,j) + settings_.g[0]);
          }
        }
    }  
    for (int j = vJBegin()+1; j < vJEnd()-1; j++)
    {
      for (int i = vIBegin()+1; i < vIEnd()-1; i++)
      {
       if (type(i,j)==0)
        {  
        G(i,j) = v(i,j) + deltat * ((computeDvDx2(i,j) + computeDvDy2(i,j)) / settings_.re - computeDuvDx(i,j) - computeDv2Dy(i,j) + settings_.g[1]);
        }
      }
      
    }
     
}

