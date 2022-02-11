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
      Array2D pvektor;
      Array2D rhsVektor;
      
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

Array2D vecAdd(Array2D A, Array2D B);
Array2D vecMultScal(double a, Array2D A);



void matrixPrint(Array2D Matrix);

public:
      CG(Discretization &discretization_);
      ~CG();
      // FieldVariable calculateP();
      virtual void calculateP();
      
      
      Array2D matMul(Array2D A, Array2D B);     
      Array2D matMulVec(Array2D A, Array2D B);     
      double matMulscal(Array2D A, Array2D B);
      void setRHSVektor();

     
};


