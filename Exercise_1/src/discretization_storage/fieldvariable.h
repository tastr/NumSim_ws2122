#pragma once

#include "array2d.h"


class FieldVariable : 
      public Array2D   
{
private:
    // Array2D data_; //wegen der Vererbeung wird das nicht gebraucht
    double maximum;
    double absmaximum;
public:
    FieldVariable(std::array<int,2> size);
    //FieldVariable(float d);
    ~FieldVariable();
    // std::array<int,2> size() const; //nicht notwendig, da in parent Klasse
    void print();
    double max() ;
    // float abs(int i, int j) ; //in parent classe definiert
    double absmax();

    // interpolation function used by output writer
    double interpolateAt(double x, double y); //TODO

    int indexconvert(int i, int j) const;    
    
    // //! access the value at coordinate (i,j), declared not const, i.e. the value can be changed
    // double &operator()(int i, int j);

    // //! get the value at coordinate (i,j), declared const, i.e. it is not possible to change the value
    // double operator()(int i, int j) const;

};


