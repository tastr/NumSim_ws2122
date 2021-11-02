#pragma once

#include "array2d.h"


class FieldVariable : 
    public Array2D
{
private:
    Array2D data_;
public:
    FieldVariable(std::array<int,2> size);
    ~FieldVariable();
};


