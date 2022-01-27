#include "centraldifferences.h"

CentralDifferences::CentralDifferences(Settings settings):Discretization(settings)
  { 
  }

 CentralDifferences::~CentralDifferences()
  {
  }

//functions for numerical derivation using the central differences scheme 
double CentralDifferences::computeDu2Dx(int i, int j) const    
 {double uiphj=(u(i,j)  + u(i+1,j)) ;
  double uimhj=u(i-1,j) + u(i,j);
   return ( uiphj*uiphj - uimhj*uimhj)/(4*delta_x); 
 } 

 double CentralDifferences::computeDv2Dy(int i, int j) const
 {double vijph = (v(i,j) + v(i,j+1));
  double vijmh = (v(i,j-1) + v(i,j));
   return ( vijph*vijph - vijmh*vijmh )/(4*delta_y);
}

 double CentralDifferences::computeDuvDx(int i, int j) const
 {
    // double viphj = (v(i,j) + v(i,j+1)) * (u(i,j)  + u(i+1,j)) ;
    // double vimhj  = (v(i-1,j) + v(i-1,j+1)) * (u(i-1,j)+ u(i,j));
    // double result_alt = (viphj - vimhj )/(delta_x*4);

    double v_iphj = (v(i+1,j)+v(i,j))/2;
    double u_ijph = (u(i,j+1)+u(i,j))/2;
    double v_imhj = (v(i,j)+v(i-1,j))/2;
    double u_imhjph= (u(i-1,j+1)+u(i-1,j))/2;
    double result_neu=((v_iphj*u_ijph) - (v_imhj*u_imhjph))/delta_x;

    return result_neu;
 }

 double CentralDifferences::computeDuvDy(int i, int j) const
 {double viphj   = (v(i,j)+v(i+1,j)) * (u(i,j)+u(i,j+1));
  double viphjm1 = (v(i,j-1)+v(i+1,j-1)) * (u(i,j-1)+u(i,j));
   return (viphj-viphjm1)/(4*delta_y);
 }


// calculates the values of auxilliarz variables F and G using the specific numerical
// derivative functions of the class
void CentralDifferences::calculation() // copied it over from donorCell
{
    for (int j = uJBegin()+1; j < uJEnd()-1; j++)
    {
        for (int i = uIBegin()+1; i < uIEnd()-1; i++)
        { 
        //  if (type(i,j)==0)
          //{  
            F(i,j) = u(i,j) + deltat * ((computeDuDx2(i,j) + computeDuDy2(i,j)) / settings_.re - computeDuvDy(i,j) - computeDu2Dx(i,j) + settings_.g[0]);
          //}//else
          //{
          //setObstacleVelocity(i,j) only nedded in one loop alternatively use loop in the function
          //}
    }  
    for (int j = vJBegin()+1; j < vJEnd()-1; j++)
    {
      for (int i = vIBegin()+1; i < vIEnd()-1; i++)
      {
     //   if (type(i,j)==0)
      //  {  
        G(i,j) = v(i,j) + deltat * ((computeDvDx2(i,j) + computeDvDy2(i,j)) / settings_.re - computeDuvDx(i,j) - computeDv2Dy(i,j) + settings_.g[1]);
      //  }
      }
      
    }
     
}

