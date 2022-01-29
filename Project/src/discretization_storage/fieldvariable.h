#pragma once

#include "array2d.h"
#include "settings.h"

class FieldVariable : public Array2D
{
private:
    double maximum;
    double absmaximum;
    std::array<double, 2> offset_;
    Settings settings_;

public:
    FieldVariable(std::array<int, 2> size, Settings settings, std::array<double, 2> offset);
    ~FieldVariable();
    void print();
    double max();
    double absmax();

    // interpolation function used by output writer
    double interpolateAt(double x, double y);

    // returns the std array index of a array2D variable
    int indexconvert(int i, int j) const;
};
