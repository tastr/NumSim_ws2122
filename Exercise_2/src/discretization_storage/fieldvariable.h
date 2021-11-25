#pragma once

#include "array2d.h"
#include "settings.h"

class FieldVariable : 
      public Array2D   
{
private:
    double maximum;
    double absmaximum;
    std::array<double,2> offset_;
    // Settings settings_;
    double dx_;
    double dy_;
public:
    FieldVariable(std::array<int,2> size, std::array<double,2> offset, std::array<double,2> mesh_width);
       ~FieldVariable();
    void print();
    double max() ;
    double absmax();

    // interpolation function used by output writer
    double interpolateAt(double x, double y); 

    //returns the std array index of a array2D variable 
    int indexconvert(int i, int j) const;    
    
    
};


