#pragma once

#include <vector>
#include <array>
#include "discretization_storage/discretization.h"
#include "settings.h"

class PressureSolver
{
protected:
     Discretization &discretization_;
     Settings settings_;

public:
     PressureSolver(Discretization &discretization);
     virtual ~PressureSolver();
     double abs_(double number);
     void setPressureBoundaries();
     void calculateRHS();
     virtual void calculateP();

     double residuum();
};
