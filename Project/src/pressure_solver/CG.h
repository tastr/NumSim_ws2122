#pragma once

#include <vector>
#include <array>
#include "pressuresolver.h"
#include "discretization_storage/array2d.h"

class CG : public PressureSolver
{
protected:
      //   double omega; //Omega festgelegt durch parameter
      Array2D Matrix;
      void setMatrix();
 
      bool isEdge(int i, int j) ;
      bool isBorder(int i, int j) ;

      void setEdge(int i, int j, int n);
      void setBorder(int i, int j, int n);
      void setInnerValue(int i, int j,int n);


void setValueTop(int i, int j, int n);
void setValueLeft(int i, int j, int n);
void setValueBottom(int i, int j, int n);
void setValueRight(int i, int j, int n);

void setUpperLeftEdge(int i, int j, int n);
void setUpperRightEdge(int i, int j, int n);
void setLowerLeftEdge(int i, int j, int n);
void setLowerRightEdge(int i, int j, int n);

public:
      CG(Discretization &discretization_);
      ~CG();
      // FieldVariable calculateP();
      virtual void calculateP();
     

     
};
