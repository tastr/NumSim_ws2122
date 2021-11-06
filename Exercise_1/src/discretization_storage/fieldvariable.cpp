#include "fieldvariable.h"
#include <iostream>


//FieldVariable::FieldVariable(std::array<int,2> size) : 
FieldVariable::FieldVariable(std::array<int,2> size)
: Array2D(size)
{
    maximum=0;
    absmaximum=0;
}

FieldVariable::~FieldVariable()
{
}

// std::array<int,2> FieldVariable::size() const
// {
//   return size();
// }

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


// double &FieldVariable::operator()(int i, int j)
// {
//   return data_(i,j);
// }

// double FieldVariable::operator()(int i, int j) const
// {
//   return data_(i,j);
// }

//Erst ueber j dann ueber i zu gehen ist effizienyter, aber fuer die Ausgabe
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

double FieldVariable::interpolateAt(double x, double y)
{
   //TODO
}