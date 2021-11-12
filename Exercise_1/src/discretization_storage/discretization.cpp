#include "discretization.h"
#include "staggeredgrid.h"
#include "settings.h"
#include "fieldvariable.h"

//constructor
Discretization::Discretization(Settings settings):
StaggeredGrid(settings) 
,F({settings.nCells[0]+1,settings.nCells[1]+2},settings,{0,0.5}) // is actually smaller than this, but makes handling easier
,G({settings.nCells[0]+2,settings.nCells[1]+1},settings,{0.5,0}) // is actually smaller than this, but makes handling easier
,rhs_({settings.nCells[0]+2,settings.nCells[1]+2},settings,{0.5,0.5}) // is actually smaller than this, but makes handling easier
{

}


// setBorderVelocity(settings.dirichletBcTop,settings.dirichletBcLeft,settings.dirichletBcRight,settings.dirichletBcBottom)
void Discretization::updateDeltaT()
{
  deltat = min2(min3(min2(delta_x,delta_y)*min2(delta_x,delta_y)*settings_.re/4,delta_x/velocity_X.absmax(),delta_y/velocity_Y.absmax())*settings_.tau, settings_.maximumDt);
}

//destructor
Discretization::~Discretization()
{
}

 double Discretization::computeDu2Dx(int i, int j) const {return 0;}   

 double Discretization::computeDv2Dy(int i, int j) const {return 0;}
 
 double Discretization::computeDuvDx(int i, int j) const {return 0;}
 
 double Discretization::computeDuvDy(int i, int j) const {return 0;}
 



// second order numerical derivation terms
double Discretization::computeDuDx2(int i, int j) const
{
 return (u(i+1,j) - 2 * u(i,j) + u(i-1,j))/(delta_x*delta_x);
}
double Discretization::computeDvDy2(int i, int j) const
{
 return (v(i,j+1) - 2 * v(i,j) + v(i,j-1))/ (delta_y*delta_y);
}

double Discretization::computeDuDy2(int i, int j) const
{
 return (u(i,j+1) - 2 * u(i,j) + u(i,j-1))/(delta_y*delta_y);
}
double Discretization::computeDvDx2(int i, int j) const
{
 return (v(i+1,j) - 2 * v(i,j) + v(i-1,j))/ (delta_x* delta_x);       
}


// first order numerical derivation terms
double Discretization::computeDuDx(int i, int j) const
{ 
  return (u(i,j)-u(i-1,j))/delta_x;
} 
double Discretization::computeDvDy(int i, int j) const
{ 
  return (v(i,j)-v(i,j-1))/ delta_y;
} 

double Discretization::computeDpDx(int i, int j) const
{
  return (p(i+1, j)-p(i,j))/delta_x;
}

double Discretization::computeDpDy(int i, int j) const
{
  return (p(i, j+1)-p(i,j))/delta_y;
}


void Discretization::updateVelocity()
{   // Does not go over Boundary
    for (int j = uJBegin()+1; j < uJEnd()-1; j++)
    {
        for (int i = uIBegin()+1; i < uIEnd()-1; i++)
        {
            velocity_X(i, j)=F(i, j)-deltat*computeDpDx(i, j);
        }
        
    }
    for (int j = vJBegin()+1; j < vJEnd()-1; j++)
    {
        for (int i = vIBegin()+1; i < vIEnd()-1; i++)
        {
            velocity_Y(i, j)=G(i, j)-deltat*computeDpDy(i, j);
        }
        
    }
    
}

void Discretization::updateBoundaryFG()
{
    for (int j = uJBegin(); j < uJEnd(); j++)
   {
        F(0,j)=u(0,j);
   }

   for (int i = vIBegin(); i < vIEnd()-1; i++)
   {
        G(i,0)=v(i,0);
   }
}
        

double Discretization::min2(double value1, double value2) const
{
    if (value1<=value2)
    {
        return value1;
    }else 
    {
        return value2;
    }
}

double Discretization::min3(double value1, double value2, double value3) const
{
    if (value1<=value2 && value1<=value3)
    {
        return value1;
    }else if (value2<=value1 && value2<=value3)
    {
        return value2;
    } else
    {
        return value3;
    }
}

double Discretization::f(int i, int j) const
{
    return F(i,j);
}
double Discretization::g(int i, int j) const
{
    return G(i,j);
}
double Discretization::rhs(int i, int j) const//const statment removed because of error
{
    return rhs_(i, j);
}
double Discretization::getOmega() const
{
   return settings_.omega; 
} 

double Discretization::getDeltaT() const
{
   return deltat;
}


void Discretization::calculation()
{
     //this should never be called, as it is a virtual function
    assert(false);

   
}

void Discretization::setRHS(int i, int j, double value)
{
    rhs_(i, j)=value;
}



