#pragma once

#include "array2d.h"


class FieldVariable : 
      public Array2D   
{
private:
    double maximum;
    double absmaximum;
public:
    FieldVariable(std::array<int,2> size);
       ~FieldVariable();
    void print();
    double max() ;
    double absmax();

    // interpolation function used by output writer
    double interpolateAt(double x, double y); 

    //returns the std array index of a array2D variable 
    int indexconvert(int i, int j) const;    
    
    
};


