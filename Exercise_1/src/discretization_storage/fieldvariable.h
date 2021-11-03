#pragma once

#include "array2d.h"


class FieldVariable : 
      public Array2D   
{
private:
    Array2D data_;
    float maximum;
    float absmaximum;
public:
    FieldVariable(std::array<int,2> size);
    //FieldVariable(float d);
    ~FieldVariable();
    std::array<int,2> size() const;
    void print();
    float max() ;
    
    //! access the value at coordinate (i,j), declared not const, i.e. the value can be changed
    double &operator()(int i, int j);

    //! get the value at coordinate (i,j), declared const, i.e. it is not possible to change the value
    double operator()(int i, int j) const;

};


