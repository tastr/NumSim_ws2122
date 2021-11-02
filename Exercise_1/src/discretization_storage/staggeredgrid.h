#pragma once

#include "fieldvariable.h"

class StaggeredGrid
{
private:
    FieldVariable pressure;
    FieldVariable velocity_X;
    FieldVariable velocity_Y;
    // maybe include A,B (and RHS)
    FieldVariable F;
    FieldVariable G;
public:
    StaggeredGrid(/* args */);
    ~StaggeredGrid();
};


