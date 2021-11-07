#pragma once

#include "discretization.h"

class DonorCell :
    public Discretization
{ 
private:
   
public:
       DonorCell(Settings settings);
       ~DonorCell();
    double computeDu2Dx(int i, int j) const;
    double computeDv2Dy(int i, int j) const;
    double computeDuvDx(int i, int j) const;
    double computeDuvDy(int i, int j) const;
};


