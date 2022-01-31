#pragma once

#include <iostream>
#include <array>
#include "discretization_storage/IntArray2d.h"

/** All settings that parametrize a simulation run.
 */
class Geometry
{
private:
  int rownumber;
  int columnnumber;
 
  int getRowNumber(std::string filename);
  int getColumnNumber(std::string filename);
  int getObstacleCount(std::string filename);
  void writeMatrix(std::string filename);
  void createGeometry();
  IntArray2D matrix;
  IntArray2D geometry;
  IntArray2D obstacleCellsIndices;
  IntArray2D fluidCellsIndices;

public:
   Geometry(std::string filename);
  ~Geometry();

    
void  printMatrix();
void printGeometry();
std::array<int, 2> getFluidCellsIndices(int i) const;
int getLengthFluidCellsIndices(int i) const;
std::array<int, 2> getObstacleCellsIndices(int i) const;
int getLengthObstacleCellsIndices(int i) const;



};
