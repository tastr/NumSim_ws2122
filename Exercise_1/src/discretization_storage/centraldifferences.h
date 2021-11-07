#pragma once

#include "discretization.h"

class CentralDifferences :
    public Discretization
{
 private:
     /* data */
 public:
    CentralDifferences(Settings settings);
     ~CentralDifferences();
    double computeDu2Dx(int i, int j) const;
    double computeDv2Dy(int i, int j) const;
    double computeDuvDx(int i, int j) const;
    double computeDuvDy(int i, int j) const;
    
    virtual void calculation();
};


