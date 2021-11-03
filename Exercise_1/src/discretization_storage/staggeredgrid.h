#pragma once

#include "fieldvariable.h"

class StaggeredGrid 
{
private:
    FieldVariable pressure, velocity_X, velocity_Y,F, G;
public:
    StaggeredGrid(std::array<int,2> size);
    ~StaggeredGrid();
     void setBorderVelocity(std::array<double,2> top,std::array<double,2> left,std::array<double,2>  right,std::array<double,2> bottom); 
     void print(std::string str);
};


