#include "fieldvariable.h"
#include <iostream>
#include <cassert>


//FieldVariable::FieldVariable(std::array<int,2> size) : 
FieldVariable::FieldVariable(std::array<int,2> size,Settings settings,std::array<double,2> offset)
: Array2D(size), offset_(offset),settings_(settings)
{
    maximum=0;
    absmaximum=0;
}

FieldVariable::~FieldVariable()
{
}

void FieldVariable::print()
{
   
   
   for (int i = 0; i < size_[0]; i++)
   {
       for (int j = 0; j < size_[1]; j++)
       {
            std::cout<< " " << (*this)(i,j) << " ";
       }
       std::cout << std::endl;
   }
}

double FieldVariable::max() 
{  
 for (int i = 0; i < size_[0]; i++)
 {
    for (int j = 0; j <  size_[1]; j++)
    {
        if (maximum<(*this)(i,j))
        {
            maximum=(*this)(i,j);  
        }
    }   
 }
 return maximum;
}

double FieldVariable::absmax() 
{  
 for (int i = 0; i < size_[0]; i++)
 {
    for (int j = 0; j <  size_[1]; j++)
    {
        if (absmaximum<abs(i,j))
        {
            absmaximum=abs(i,j);
        }        
    }   
 }
 return absmaximum;
}

int FieldVariable::indexconvert(int i, int j) const
{
      // assert that indices are in range
     assert(0 <= i && i < size_[0]);
     assert(0 <= j && j < size_[1]);
     assert(j*size_[0] + i < (int)data_.size());
    return  j*size_[0] + i;
}


double FieldVariable::interpolateAt(double x, double y)
{ int i_right=1;
  int j_top=1;
  double propfact_x;
  double propfact_y;
  double top_average; 
  double bottom_average;
  double averaged=0;
  
  double dx=settings_.physicalSize[0] / (1.0*settings_.nCells[0]);
  double dy=settings_.physicalSize[1] / (1.0*settings_.nCells[1]);

  x=x+ offset_[0]*dx;
  y=y+ offset_[1]*dy; 
 

    while (x>(i_right*dx))
    {
        i_right=++i_right;
       // std::cout << i_right <<std::endl;
    }
    while (y>(dy*j_top))
    { 
        j_top=++j_top;     
    }
 
    // propfact_x =x*settings_.nCells[0]-i_right+1;
    propfact_x = (x-(i_right-1)*dx)/dx;
    // propfact_y =y*settings_.nCells[1]-j_top+1;
    propfact_y = (y-(j_top-1)*dy)/dy;
 
//    std::cout << propfact_x <<"  " <<  propfact_y <<std::endl;
  bottom_average = (1-propfact_x) * (*this)(i_right-1,j_top-1) + (propfact_x) * (*this)(i_right,j_top-1);
  top_average    = (1-propfact_x) * (*this)(i_right-1,j_top)   + (propfact_x) * (*this)(i_right,j_top);
  averaged       = (1-propfact_y) * bottom_average + propfact_y * top_average;
  //std::cout << i_right << "   " << x <<std::endl;
  //std::cout << j_top << "   " << y <<std::endl;
//   std::cout << "outupt of interplateAt: "<< averaged <<std::endl;
  return  averaged ;
//    return 1;
}