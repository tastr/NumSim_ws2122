#pragma once

#include <iostream>
#include <array>
#include "discretization_storage/IntArray2d.h"


class Geometry
{
private:
    int rownumber;
    int columnnumber;
    // Matrix with the values from the file
    IntArray2D matrix;
    // Array that contains the flags of the obsticales
    IntArray2D geometry;
    // Array with the indices of the obstacle Values
    IntArray2D obstacleCellsIndices;
    // Array with the indices of the fluid Values
    IntArray2D fluidCellsIndices;

    /*
     * returns the number of rows in the file
     */
    int getRowNumber(std::string filename);

    /*
     * returns the number of columns in the file
     */
    int getColumnNumber(std::string filename);

    /*
     * returns the number of obstacles based on the original file (for pre constructor initialisation)
     */
    int getObstacleCount(std::string filename);

    /*
     * writes the values from the file to an array
     */
    void writeMatrix(std::string filename);

    /*
     * creates the geometry based on the values of the matrix
     */
    void createGeometry();

public:
    /*
     * constructor, gets thye name of the file containming the matrix with the geometry
     */
    Geometry(std::string filename);

    /*
     * destructor
     */
    ~Geometry();

    /*
     * prints the Values whom have been read from the file onto the consol
     * Starts counting in the top left corner
     */
    void printMatrix();

    /*
     * prints the Values of the geometry
     * Starts counting in the bottom left corner
     */
    void printGeometry();

    /*
     * returns the indices of the ith fluid cell out of the vector with fluidcells
     */
    std::array<int, 2> getFluidCellsIndices(int i) const;

    /*
     * returns the length of the fluidcellvector (number of fluidcells)
     */
    int getLengthFluidCellsIndices(int i) const;

    /*
     * returns the indices of the ith obstacles cell out of the vector with fluidcells
     */
    std::array<int, 2> getObstacleCellsIndices(int i) const;

    /*
     * returns the length of the obstaclecellvector (number of obstaclecells)
     */
    int getLengthObstacleCellsIndices(int i) const;

    IntArray2D getGeometry() const;
};
