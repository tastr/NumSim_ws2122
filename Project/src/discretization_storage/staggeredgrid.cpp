#include "staggeredgrid.h"
#include "settings.h"
#include <iostream>
#include <cassert>
#include "Geometry.h"

// StaggeredGrid::StaggeredGrid(std::array<int,2> size)
StaggeredGrid::StaggeredGrid(Settings settings)
    : pressure({settings.nCells[0] + 2, settings.nCells[1] + 2},
               settings, {0.5, 0.5}),
      velocity_X({settings.nCells[0] + 1, settings.nCells[1] + 2},
                 settings, {0, 0.5}),
      velocity_Y({settings.nCells[0] + 2, settings.nCells[1] + 1},
                 settings, {0.5, 0}),
      type({settings.nCells[0] + 2, settings.nCells[1] + 2}, settings, {0, 0}),
      settings_(settings)
//: pressure({settings.nCells[0]+2,settings.nCells[1]+2})
//,velocity_X({settings.nCells[0]+2,settings.nCells[1]+2})
//,velocity_Y({settings.nCells[0]+2,settings.nCells[1]+2})
//,settings_(settings)

{
    setSize_(settings.nCells);
    delta_x = settings_.physicalSize[0] / (1.0 * settings_.nCells[0]);
    delta_y = settings_.physicalSize[1] / (1.0 * settings_.nCells[1]);
    epsilon = settings.epsilon;

    uOldLeft.resize(size_[1], 0.0);
    uOldRight.resize(size_[1], 0.0);
    vOldTop.resize(size_[0], 0.0);
    vOldBottom.resize(size_[0], 0.0);
}

void StaggeredGrid::setSize_(std::array<int, 2> nCells)
{
    size_ = {nCells[0] + 2, nCells[1] + 2};
}

std::array<int, 2> StaggeredGrid::getSize() const
{
    return size_;
}

int StaggeredGrid::getMaxIteration() const
{
    return settings_.maximumNumberOfIterations;
}

double StaggeredGrid::getURrhs() const
{
    return settings_.underrelaxationRHS;
}

// Sets the Dirichlet border conditions for the velocity variables
void StaggeredGrid::setBorderVelocity(std::array<double, 2> top, std::array<double, 2> left, std::array<double, 2> right, std::array<double, 2> bottom)
{
    if (settings_.outflowLeft == true)
    {
        for (int j = 0; j < uJEnd(); j++)
        {
            velocity_X(uIBegin(), j) = u(uIBegin() + 1, j);
        }
        for (int j = vJBegin(); j < vJEnd(); j++)
        {
            velocity_Y(vIBegin(), j) = velocity_Y(vIBegin() + 1, j);
        }
    } 
    else
    {
        for (int j = 0; j < uJEnd(); j++)
        {
            velocity_X(uIBegin(), j) = left[0];
        }
        for (int j = vJBegin(); j < vJEnd(); j++)
        {
            velocity_Y(vIBegin(), j) = 2 * left[1] - velocity_Y(vIBegin() + 1, j);
        }
    }

    if (settings_.outflowRight == true)
    {
        for (int j = 0; j < uJEnd(); j++)
        {
            velocity_X(uIEnd() - 1, j) = u(uIEnd() - 2, j);
        }
        for (int j = vJBegin(); j < vJEnd(); j++)
        {
            velocity_Y(vIEnd() - 1, j) = velocity_Y(vIEnd() - 2, j);
        }
    } 
    else
    {
        for (int j = 0; j < uJEnd(); j++)
        {
            velocity_X(uIEnd() - 1, j) = right[0];
        }
        for (int j = vJBegin(); j < vJEnd(); j++)
        {
            velocity_Y(vIEnd() - 1, j) = 2 * right[1] - velocity_Y(vIEnd() - 2, j);
        }
    }

    if (settings_.outflowBottom == true)
    {
        for (int i = 1; i < uIEnd() - 1; i++) // i starts at 1 and goes to i_u_max-1 so that the wall is the BC in corners
        {
            velocity_X(i, uJBegin()) = velocity_X(i, vJBegin() + 1);
        }
        for (int i = vIBegin() + 1; i < vIEnd() - 1; i++) // i starts at 1 and goes to i_u_max-1 so that the wall is the BC in corners
        {
            velocity_Y(i, vJBegin()) = velocity_Y(i, vJBegin() + 1);
        }
    } 
    else
    {
        for (int i = 1; i < uIEnd() - 1; i++) // i starts at 1 and goes to i_u_max-1 so that the wall is the BC in corners
        {
            velocity_X(i, uJBegin()) = 2 * bottom[0] - velocity_X(i, vJBegin() + 1);
        }
        for (int i = vIBegin() + 1; i < vIEnd() - 1; i++) // i starts at 1 and goes to i_u_max-1 so that the wall is the BC in corners
        {
            velocity_Y(i, vJBegin()) = bottom[1];
        }
    }

    if (settings_.outflowTop == true)
    {
        for (int i = 1; i < uIEnd() - 1; i++) // i starts at 1 and goes to i_u_max-1 so that the wall is the BC in corners
        {
            velocity_X(i, uJEnd() - 1) = velocity_X(i, uJEnd() - 2);
        }
        for (int i = vIBegin() + 1; i < vIEnd() - 1; i++) // i starts at 1 and goes to i_u_max-1 so that the wall is the BC in corners
        {
            velocity_Y(i, vJEnd() - 1) = velocity_Y(i, vJEnd() - 2);
        }
    } 
    else
    {
        for (int i = 1; i < uIEnd() - 1; i++) // i starts at 1 and goes to i_u_max-1 so that the wall is the BC in corners
        {
            velocity_X(i, uJEnd() - 1) = 2 * top[0] - velocity_X(i, uJEnd() - 2);
        }
        for (int i = vIBegin() + 1; i < vIEnd() - 1; i++) // i starts at 1 and goes to i_u_max-1 so that the wall is the BC in corners
        {
            velocity_Y(i, vJEnd() - 1) = top[1];
        }
    }
    
    /*
    for (int j = 0; j < uJEnd(); j++)
    {
        velocity_X(uIBegin(), j) = left[0];
        velocity_X(uIEnd() - 1, j) = right[0];
    }
    for (int i = 1; i < uIEnd() - 1; i++) // i starts at 1 and goes to i_u_max-1 so that the wall is the BC in corners
    {
        velocity_X(i, uJBegin()) = 2 * bottom[0] - velocity_X(i, vJBegin() + 1);
        velocity_X(i, uJEnd() - 1) = 2 * top[0] - velocity_X(i, uJEnd() - 2);
    }

    // set v velocity
    for (int j = vJBegin(); j < vJEnd(); j++)
    {
        velocity_Y(vIBegin(), j) = 2 * left[1] - velocity_Y(vIBegin() + 1, j);
        velocity_Y(vIEnd() - 1, j) = 2 * right[1] - velocity_Y(vIEnd() - 2, j);
    }
    for (int i = vIBegin() + 1; i < vIEnd() - 1; i++) // i starts at 1 and goes to i_u_max-1 so that the wall is the BC in corners
    {
        velocity_Y(i, vJBegin()) = bottom[1];
        velocity_Y(i, vJEnd() - 1) = top[1];
    }
    */
}

// updates the border terms for the pressure variable in the end of each iteration
void StaggeredGrid::updatedPressureBC()
{
    int i_max = pressure.size()[0], j_max = pressure.size()[1];
    for (int j = 0; j < j_max; j++)
    {
        pressure(0, j) = pressure(1, j);
        pressure(i_max - 1, j) = pressure(i_max - 2, j);
    }
    for (int i = 0; i < i_max; i++)
    {
        pressure(i, 0) = pressure(i, 1);
        pressure(i, j_max - 1) = pressure(i, j_max - 2);
    }
}

// prints out the values stored in the FieldVariables
void StaggeredGrid::print(std::string str)
{
    if (str == "pressure")
    {
        pressure.print();
    }
    else if (str == "velocity_X" || str == "velocity_x")
    {
        velocity_X.print();
    }
    else if (str == "velocity_Y" || str == "velocity_y")
    {
        velocity_Y.print();
    }
    else if (str == "all")
    {
        std::cout << "pressure" << std::endl;
        pressure.print();
        std::cout << "velocity_X" << std::endl;
        velocity_X.print();
        std::cout << "velocity_Y" << std::endl;
        velocity_Y.print();
    }
}

// functions return values describing the grid
std::array<double, 2> StaggeredGrid::meshWidth() const
{
    return {dx(), dy()};
}
std::array<int, 2> StaggeredGrid::nCells() const
{
    return settings_.nCells;
}
double StaggeredGrid::getepsilon() const
{
    return epsilon;
}
double StaggeredGrid::dx() const
{
    return delta_x;
}
double StaggeredGrid::dy() const
{
    return delta_y;
}

// functions to acess entire Fieldvariables of the grid
FieldVariable StaggeredGrid::p() const
{
    return pressure;
}
FieldVariable StaggeredGrid::u() const
{
    return velocity_X;
}
FieldVariable StaggeredGrid::v() const
{
    return velocity_Y;
}

// functions to retrieve singular values from the FieldVariables
double StaggeredGrid::u(int i, int j) const
{
    return velocity_X(i, j);
}
double StaggeredGrid::p(int i, int j) const
{
    return pressure(i, j);
}

double StaggeredGrid::v(int i, int j) const
{
    return velocity_Y(i, j);
}

double StaggeredGrid::getTyp(int i, int j) const
{
    return type(i, j);
}

// functions return the size values of the FieldVariables
int StaggeredGrid::uIBegin() const
{
    return 0;
    // return 1;
}
int StaggeredGrid::uIEnd() const
{
    return velocity_X.size()[0];
}
int StaggeredGrid::uJBegin() const
{
    return 0;
}
int StaggeredGrid::uJEnd() const
{
    return velocity_X.size()[1];
}
int StaggeredGrid::vIBegin() const
{
    return 0;
}
int StaggeredGrid::vIEnd() const
{
    return velocity_Y.size()[0];
}
int StaggeredGrid::vJBegin() const
{
    return 0;
    // return 1;
}
int StaggeredGrid::vJEnd() const
{
    return velocity_Y.size()[1];
}
int StaggeredGrid::pIBegin() const
{
    return 0;
}
int StaggeredGrid::pIEnd() const
{
    return pressure.size()[0];
}
int StaggeredGrid::pJBegin() const
{
    return 0;
}
int StaggeredGrid::pJEnd() const
{
    return pressure.size()[1];
}

// functions to write values into the Fieldvariables
void StaggeredGrid::setU(int i, int j, double value)
{
    velocity_X(i, j) = value;
}
void StaggeredGrid::setV(int i, int j, double value)
{
    velocity_Y(i, j) = value;
}
void StaggeredGrid::setP(int i, int j, double value)
{
    pressure(i, j) = value;
}

// functions to overwrite the FieldVariables of the grid with the Values of other Fieldvariables with maching size
void StaggeredGrid::setV(FieldVariable value)
{ // assert that indices are in range
    assert(value.size()[0] == velocity_Y.size()[0]);
    assert(value.size()[1] == velocity_Y.size()[1]);
    for (int j = 0; j < value.size()[1]; j++)
    {
        for (int i = 0; i < value.size()[0]; i++)
        {
            velocity_Y(i, j) = value(i, j);
        }
    }
}

void StaggeredGrid::setP(FieldVariable value)
{ // assert that indices are in range
    assert(value.size()[0] == pressure.size()[0]);
    assert(value.size()[1] == pressure.size()[1]);
    for (int j = 0; j < value.size()[1]; j++)
    {
        for (int i = 0; i < value.size()[0]; i++)
        {
            pressure(i, j) = value(i, j);
        }
    }
}

// function to calculate the absolute value of a double
double StaggeredGrid::abs(double number) const
{
    if (number >= 0)
    {
        return number;
    }
    else
    {
        return -number;
    }
}



void StaggeredGrid::setObstaclePressure()
{
    for (int j = 1; j < type.size()[1] - 1; j++)
    {
        for (int i = 1; i < type.size()[0] - 1; i++)
        {
            if (type(i, j) != 0 && type(i, j) != 1)
            {
                if (type(i, j) == 2)
                {
                    pressure(i, j) = p(i, j + 1);
                }
                else if (type(i, j) == 3)
                {
                    pressure(i, j) = p(i, j - 1);
                }
                else if (type(i, j) == 4)
                {
                    pressure(i, j) = p(i + 1, j);
                }
                else if (type(i, j) == 5)
                {
                    pressure(i, j) = p(i - 1, j);
                }
                else if (type(i, j) == 6)
                {
                    pressure(i, j) = 0.5 * (p(i + 1, j) + p(i, j + 1));
                }
                else if (type(i, j) == 7)
                {
                    pressure(i, j) = 0.5 * (p(i - 1, j) + p(i, j + 1));
                }
                else if (type(i, j) == 8)
                {
                    pressure(i, j) = 0.5 * (p(i - 1, j) + p(i, j - 1));
                }
                else if (type(i, j) == 9)
                {
                    pressure(i, j) = 0.5 * (p(i + 1, j) + p(i, j - 1));
                }
            }
        }
    }
}
void StaggeredGrid::setObstacle(Geometry geometrie)
{
    /*
    for (int i = 1; i < 7; i++)
    {
        type(i, 15) = 3;
        type(i, 16) = 2;
    }

    for (int i = 14; i < type.size()[0] - 1; i++)
    {
        type(i, 15) = 3;
        type(i, 16) = 2;
    }

    type(7, 15) = 9;
    type(7, 16) = 6;

    type(13, 15) = 8;
    type(13, 16) = 7;
    */

/*
 for (int i = 1; i < 7; i++)
    {
        type(i, 13) = 2;
        type(i, 18) = 3;
       
        type(i, 14) = 1;
        type(i, 15) = 1;
        type(i, 16) = 1;
        type(i, 17) = 1;
    }

    for (int i = 14; i < type.size()[0] - 1; i++)
    {
        type(i, 13) = 2;
        type(i, 18) = 3;

        type(i, 14) = 1;
        type(i, 15) = 1;
        type(i, 16) = 1;
        type(i, 17) = 1;
    }

    type(7, 13) = 6;
    type(7, 18) = 9;

    type(7, 14) = 4;
    type(7, 15) = 4;
    type(7, 16) = 4;
    type(7, 17) = 4;


    type(13, 13) = 7;
    type(13, 18) = 8;

    type(13, 14) = 5;
    type(13, 15) = 5;
    type(13, 16) = 5;
    type(13, 17) = 5;
*/
/*

 for (int i = 1; i < 8; i++)
    {
        type(i, 12) = 2;
        type(i, 19) = 3;
       
        type(i, 14) = 1;
        type(i, 15) = 1;
        type(i, 16) = 1;
        type(i, 17) = 1;
        type(i, 13) = 1;
        type(i, 18) = 1;
    }

    for (int i = 13; i < type.size()[0] - 1; i++)
    {
        type(i, 12) = 2;
        type(i, 19) = 3;

        type(i, 14) = 1;
        type(i, 15) = 1;
        type(i, 16) = 1;
        type(i, 17) = 1;
        type(i, 13) = 1;
        type(i, 18) = 1;
        
        

    }

    type(8, 19) = 6;
    type(8, 12) = 9;

    type(8, 14) = 4;
    type(8, 15) = 4;
    type(8, 16) = 4;
    type(8, 17) = 4;
    type(8, 13) = 4;
    type(8, 18) = 4;
    

    type(12, 19) = 7;
    type(12, 12) = 8;

    type(12, 14) = 5;
    type(12, 15) = 5;
    type(12, 16) = 5;
    type(12, 17) = 5;
    type(12, 13) = 5;
    type(12, 18) = 5;

*/


// 30x20 grid
/*
    for (int i = 1; i < 7; i++)
    {
        type(15, i) = 2;
        type(16, i) = 3;
    }

    for (int i = 14; i < type.size()[1] - 1; i++)
    {
        type(15, i) = 2;
        type(16, i) = 3;
    }

    type(15, 7) = 7;
    type(16, 7) = 6;

    type(15, 13) = 8;
    type(16, 13) = 9;
*/

/*
  for (int i = 1; i < 7; i++)
    {   
        type(14, i) = 5;
        type(15, i) = 1;
        type(16, i) = 1;
        type(17, i) = 4;

    }

    type(14, 7) = 7;
    type(15, 7) = 2;
    type(16, 7) = 2;
    type(17, 7) = 6;



 
type.print();

*/







printf("%d %d",geometrie.getGeometry().size()[0],geometrie.getGeometry().size()[1]);

IntArray2D geo=geometrie.getGeometry();
for (int i = 0; i<type.size()[0] ; i++)
{
    for (int j = 0; j<type.size()[1] ; j++)
    {
     type(i,j)=geo(i,j);
    }
    
}






}

// TODO setObstacleFG