#include "fieldvariable.h"
#include <iostream>


//FieldVariable::FieldVariable(std::array<int,2> size) : 
FieldVariable::FieldVariable(std::array<int,2> size)
: Array2D(size), data_(size)
{
    maximum=0;
    absmaximum=0;
}

FieldVariable::~FieldVariable()
{
}

std::array<int,2> FieldVariable::size() const
{
  return data_.size();
}

void FieldVariable::print()
{
   
   std::array<int,2> size=data_.size();
   for (int j = 0; j < size[1]; j++)
   {
       for (int i = 0; i < size[0]; i++)
       {
            std::cout<< " " << data_(i,j) << " ";
       }
       std::cout << std::endl;
   }
}


double &FieldVariable::operator()(int i, int j)
{
  return data_(i,j);
}

double FieldVariable::operator()(int i, int j) const
{
  return data_(i,j);
}

float FieldVariable::max() 
{
 for (int j = 0; j < data_.size()[1]; j++)
 {
    for (int i = 0; i <  data_.size()[0]; i++)
    {
        if (maximum<data_(i,j))
        {
            maximum=data_(i,j);
        }
    }   
 }
 return maximum;
}
