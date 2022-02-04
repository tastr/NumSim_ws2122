#include <fstream>
#include <iomanip>
#include <stdexcept>
#include "Geometry.h"

Geometry::Geometry(std::string filename)
    : matrix(IntArray2D({getColumnNumber(filename), getRowNumber(filename)})),
      geometry(IntArray2D({getColumnNumber(filename) + 2, getRowNumber(filename) + 2})),
      fluidCellsIndices(IntArray2D({2, getColumnNumber(filename) * getRowNumber(filename) - getObstacleCount(filename)})),
      obstacleCellsIndices(IntArray2D({2, getObstacleCount(filename)}))
{
    rownumber = getRowNumber(filename);
    columnnumber = getColumnNumber(filename);
    writeMatrix(filename);
    createGeometry();
}
Geometry::~Geometry()
{
}

int Geometry::getRowNumber(std::string filename)
{
    std::ifstream file(filename.c_str(), std::ios::in);
    if (!file.is_open())
    {
        std::cout << "Could not open Geometrie file \"" << filename << "\"." << std::endl;
        throw std::runtime_error("Could not open Geometrie file");
    }
    else
    {
        int rownumber = -1;

        while (!file.eof())
        {
            std::string line;
            getline(file, line);
            rownumber++;
        }
        return rownumber;
    }
    
}

int Geometry::getColumnNumber(std::string filename)
{
    std::ifstream file(filename.c_str(), std::ios::in);
    if (!file.is_open())
    {
        std::cout << "Could not open Geometrie file \"" << filename << "\"." << std::endl;
        throw std::runtime_error("Could not open Geometrie file");
    }
        else
        {
        int columnnumber = 0;
        std::string line;

        getline(file, line);
        while (line.find(",") != std::string::npos)
        {
            columnnumber++;
            line.erase(line.find_first_of(","), 1);
        }

        return columnnumber;
    }
}

void Geometry::writeMatrix(std::string filename)
{
    std::ifstream file(filename.c_str(), std::ios::in);
    if (!file.is_open())
    {
        std::cout << "Could not open Geometrie file \"" << filename << "\"." << std::endl;
        throw std::runtime_error("Could not open Geometrie file");
    }
    else
    {
        int i = 0;
        int j = 0;
        
        std::string line;
        std::string parameter = "";

        while (!file.eof())
        {
            getline(file, line);
            // TODO Solution for the Problem, that , at the end of line is required
            while (line.find(",") != std::string::npos)
            {

                parameter = line.substr(0, line.find_first_of(","));
                // matrix(i, j) = atof(parameter.c_str());
                matrix(i, matrix.size()[1] - 1 - j) = atoi(parameter.c_str());
                line.erase(0, line.find_first_of(",") + 1);
                i++;
            }
            i = 0;
            j++;
        }
    }
}

void Geometry::printMatrix()
{
    for (int j = 0; j < matrix.size()[1]; j++)
    {
        for (int i = 0; i < matrix.size()[0]; i++)
        {
            printf("%d ", matrix(i, matrix.size()[1] - 1 - j));
        }
        printf("\n");
    }
    printf("\n");
}

void Geometry::printGeometry()
{
    for (int j = 0; j < geometry.size()[1]; j++)
    {
        for (int i = 0; i < geometry.size()[0]; i++)
        {
            printf("%d ", geometry(i, rownumber + 1 - j));
        }
        printf("\n");
    }
    printf("\n");
}

int Geometry::getObstacleCount(std::string filename)
{
    int obstaclecount = 0;
    std::ifstream file(filename.c_str(), std::ios::in);
    std::string line;
    std::string parameter = "";

    while (!file.eof())
    {
        getline(file, line);
        // TODO Solution for the Problem, that , at the end of line is required
        while (line.find(",") != std::string::npos)
        {

            parameter = line.substr(0, line.find_first_of(","));
            if (atoi(parameter.c_str()) == 1)
            {
                obstaclecount++;
            }
            line.erase(0, line.find_first_of(",") + 1);
        }
    }
    return obstaclecount;
}

void Geometry::createGeometry()
{
    int fluidindex = 0;
    int obstacleindex = 0;
    int sumisfluidcell = 0;

    for (int j = 0; j < matrix.size()[1]; j++)
    {
        for (int i = 0; i < matrix.size()[0]; i++)
        {
            geometry(i + 1, j + 1) = matrix(i, j);
        }
    }
    for (int j = 0; j < geometry.size()[1]; j++)
    {
        geometry(0, j) = 1;
        geometry(geometry.size()[0] - 1, j) = 1;
    }
    for (int i = 0; i < geometry.size()[0]; i++)
    {
        geometry(i, 0) = 1;
        geometry(i, geometry.size()[1] - 1) = 1;
    }

    for (int j = 1; j < geometry.size()[1] - 1; j++)
    {
        for (int i = 1; i < geometry.size()[0] - 1; i++)
        { // Clockwise Top, right, Bottom, left
            if (geometry(i, j) == 1)
            {
                int isfluidcell[4] = {geometry(i, j + 1) == 0, geometry(i + 1, j) == 0, geometry(i, j - 1) == 0, geometry(i - 1, j) == 0};
                sumisfluidcell = 0;
                for (int k = 0; k < 4; k++)
                {
                    sumisfluidcell = sumisfluidcell + isfluidcell[k];
                }
                if (sumisfluidcell > 2 || (isfluidcell[0] && isfluidcell[2]) || (isfluidcell[1] && isfluidcell[3]))
                {
                    printf("%d %d %d %d \n", isfluidcell[0], isfluidcell[1], isfluidcell[2], isfluidcell[3]);
                    printf("Koordinate i=%d und j=%d  geometry=%d \n", j, i, geometry(i, j));
                    printf("\n");
                    printMatrix();
                    throw std::invalid_argument("geometry verletzt zwei-Zellen-Regel");
                }

                obstacleCellsIndices(0, obstacleindex) = i;
                obstacleCellsIndices(1, obstacleindex) = j;
                obstacleindex++;

                if (sumisfluidcell == 1)
                {
                    if (isfluidcell[0])
                    {
                        geometry(i, j) = 2;
                    }
                    else if (isfluidcell[2])
                    {
                        geometry(i, j) = 3;
                    }
                    else if (isfluidcell[1])
                    {
                        geometry(i, j) = 4;
                    }
                    else if (isfluidcell[3])
                    {
                        geometry(i, j) = 5;
                    }
                }
                else
                {
                    if (isfluidcell[0] && isfluidcell[1])
                    {
                        geometry(i, j) = 6;
                    }
                    else if (isfluidcell[0] && isfluidcell[3])
                    {
                        geometry(i, j) = 7;
                    }
                    else if (isfluidcell[2] && isfluidcell[3])
                    {
                        geometry(i, j) = 8;
                    }
                    else if (isfluidcell[1] && isfluidcell[2])
                    {
                        geometry(i, j) = 9;
                    }
                }
            }
            else
            {
                fluidCellsIndices(0, fluidindex) = i;
                fluidCellsIndices(1, fluidindex) = j;
                fluidindex++;
            }
        }
    }
}

std::array<int, 2> Geometry::getFluidCellsIndices(int i) const
{
    return {fluidCellsIndices(0, i), fluidCellsIndices(1, i)};
}
int Geometry::getLengthFluidCellsIndices(int i) const
{
    return fluidCellsIndices.size()[1];
}
std::array<int, 2> Geometry::getObstacleCellsIndices(int i) const
{
    return {obstacleCellsIndices(0, i), obstacleCellsIndices(1, i)};
}
int Geometry::getLengthObstacleCellsIndices(int i) const
{
    return obstacleCellsIndices.size()[1];
}

IntArray2D Geometry::getGeometry() const
{
    return geometry;
}

std::array<int, 2> Geometry::getNumberOfCells() const
{
    return {columnnumber,rownumber};
}



