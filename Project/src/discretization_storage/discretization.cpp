#include "discretization.h"
#include "staggeredgrid.h"
#include "settings.h"
#include "fieldvariable.h"

// constructor
Discretization::Discretization(Settings settings) : StaggeredGrid(settings), F({settings.nCells[0] + 1, settings.nCells[1] + 2}, settings, {0, 0.5}) // is actually smaller than this, but makes handling easier
                                                    ,
                                                    G({settings.nCells[0] + 2, settings.nCells[1] + 1}, settings, {0.5, 0}) // is actually smaller than this, but makes handling easier
                                                    ,
                                                    rhs_({settings.nCells[0] + 2, settings.nCells[1] + 2}, settings, {0.5, 0.5}) // is actually smaller than this, but makes handling easier
{
}

// setBorderVelocity(settings.dirichletBcTop,settings.dirichletBcLeft,settings.dirichletBcRight,settings.dirichletBcBottom)
void Discretization::updateDeltaT()
{
    double time_limit_diffusion = (delta_x * delta_x * delta_y * delta_y) / ((delta_x * delta_x) + (delta_y * delta_y)) * settings_.re / 2;
    double time_limit = min3(time_limit_diffusion, delta_x / velocity_X.absmax(), delta_y / velocity_Y.absmax()) * settings_.tau;
    deltat = min2(time_limit, settings_.maximumDt);
}

// destructor
Discretization::~Discretization()
{
}

double Discretization::computeDu2Dx(int i, int j) const { return 0; }

double Discretization::computeDv2Dy(int i, int j) const { return 0; }

double Discretization::computeDuvDx(int i, int j) const { return 0; }

double Discretization::computeDuvDy(int i, int j) const { return 0; }

// second order numerical derivation terms
double Discretization::computeDuDx2(int i, int j) const
{
    return (u(i + 1, j) - 2 * u(i, j) + u(i - 1, j)) / (delta_x * delta_x);
}
double Discretization::computeDvDy2(int i, int j) const
{
    return (v(i, j + 1) - 2 * v(i, j) + v(i, j - 1)) / (delta_y * delta_y);
}

double Discretization::computeDuDy2(int i, int j) const
{
    return (u(i, j + 1) - 2 * u(i, j) + u(i, j - 1)) / (delta_y * delta_y);
}
double Discretization::computeDvDx2(int i, int j) const
{
    return (v(i + 1, j) - 2 * v(i, j) + v(i - 1, j)) / (delta_x * delta_x);
}

// first order numerical derivation terms
double Discretization::computeDuDx(int i, int j) const
{
    return (u(i, j) - u(i - 1, j)) / delta_x;
}
double Discretization::computeDvDy(int i, int j) const
{
    return (v(i, j) - v(i, j - 1)) / delta_y;
}

double Discretization::computeDpDx(int i, int j) const
{
    return (p(i + 1, j) - p(i, j)) / delta_x;
}

double Discretization::computeDpDy(int i, int j) const
{
    return (p(i, j + 1) - p(i, j)) / delta_y;
}

void Discretization::updateVelocity()
{ // Does not go over Boundary
    for (int j = uJBegin() + 1; j < uJEnd() - 1; j++)
    {
        for (int i = uIBegin() + 1; i < uIEnd() - 1; i++)
        {
            if (type(i, j) == 0)
            {
                velocity_X(i, j) = (1 - settings_.underrelaxationUV) * velocity_X(i, j) + settings_.underrelaxationUV * (F(i, j) - deltat * computeDpDx(i, j));
            }
        }
    }
    for (int j = vJBegin() + 1; j < vJEnd() - 1; j++)
    {
        for (int i = vIBegin() + 1; i < vIEnd() - 1; i++)
        {
            if (type(i, j) == 0)
            {
                velocity_Y(i, j) = (1 - settings_.underrelaxationUV) * velocity_Y(i, j) + settings_.underrelaxationUV * (G(i, j) - deltat * computeDpDy(i, j));
            }
        }
    }
}

void Discretization::updateBoundaryFG()
{
    if (settings_.outflowLeft == true)
    {
        for (int j = uJBegin(); j < uJEnd(); j++)
        {
            F(0, j) = 2 * u(0, j) - uOldLeft[j];
            uOldLeft[j] = u(0, j); // stroring current velocity to use in next time step
        }
    } else
    {
        for (int j = uJBegin(); j < uJEnd(); j++)
        {
            F(0, j) = u(0, j);
        }
    }

    if (settings_.outflowRight == true)
    {
        for (int j = uJBegin(); j < uJEnd(); j++)
        {
            F(uIEnd() - 1, j) = 2 * u(uIEnd() - 1, j) - uOldRight[j];
            uOldRight[j] = u(uIEnd() - 1, j);
        }
    } else
    {
        for (int j = uJBegin(); j < uJEnd(); j++)
        {
            F(uIEnd() - 1, j) = u(uIEnd() - 1, j);
        }
    }

    if (settings_.outflowBottom == true)
    {
        for (int i = vIBegin(); i < vIEnd(); i++)
        {
            G(i, vJBegin()) = 2 * v(i, vJBegin()) - vOldBottom[i];
            vOldBottom[i] = v(i, vJBegin());
        }
    } else
    {
        for (int i = vIBegin(); i < vIEnd(); i++)
        {
            G(i, vJBegin()) = v(i, vJBegin());
        }
    }

    if (settings_.outflowTop == true)
    {
        for (int i = vIBegin(); i < vIEnd(); i++)
        {
            G(i, vJEnd() - 1) = 2 * v(i, vJEnd() - 1) - vOldTop[i];
            vOldTop[i]  = v(i, vJEnd() - 1);
        }
    } else
    {
        for (int i = vIBegin(); i < vIEnd(); i++)
        {
            G(i, vJEnd() - 1) = v(i, vJEnd() - 1);
        }
    }
}

double Discretization::min2(double value1, double value2) const
{
    if (value1 <= value2)
    {
        return value1;
    }
    else
    {
        return value2;
    }
}

double Discretization::min3(double value1, double value2, double value3) const
{
    if (value1 <= value2 && value1 <= value3)
    {
        return value1;
    }
    else if (value2 <= value1 && value2 <= value3)
    {
        return value2;
    }
    else
    {
        return value3;
    }
}

double Discretization::f(int i, int j) const
{
    return F(i, j);
}
double Discretization::g(int i, int j) const
{
    return G(i, j);
}
double Discretization::rhs(int i, int j) const // const statment removed because of error
{
    return rhs_(i, j);
}
double Discretization::getOmega() const
{
    return settings_.omega;
}

double Discretization::getDeltaT() const
{
    return deltat;
}

void Discretization::calculation()
{
    // this should never be called, as it is a virtual function
    assert(false);
}

void Discretization::setRHS(int i, int j, double value)
{
    rhs_(i, j) = value;
}

// needs to be divided in loops for u and v, alternatively, no 2...9 obstacles on the bounderies
void Discretization::setObstacleVelocityFG()
{
    for (int j = 1; j < type.size()[1] - 1; j++)
    {
        for (int i = 1; i < type.size()[0] - 1; i++)
        {
            if (type(i, j) != 0 && type(i, j) != 1)
            {
                if (type(i, j) == 2)
                {
                    velocity_Y(i, j) = 0;
                    velocity_X(i, j) = -u(i, j + 1);
                    velocity_X(i - 1, j) = -u(i - 1, j + 1);
                    G(i, j) = 0;
                }
                else if (type(i, j) == 3)
                {
                    velocity_Y(i, j - 1) = 0;
                    velocity_X(i, j) = -u(i, j - 1);
                    velocity_X(i - 1, j) = -u(i - 1, j - 1);
                    G(i, j - 1) = 0;
                }
                else if (type(i, j) == 4)
                {
                    velocity_Y(i, j) = -v(i + 1, j);
                    velocity_X(i, j) = 0;
                    velocity_Y(i, j - 1) = -v(i + 1, j - 1);
                    F(i, j) = 0;
                }
                else if (type(i, j) == 5)
                {
                    velocity_Y(i, j) = -v(i - 1, j);
                    velocity_X(i - 1, j) = 0;
                    velocity_Y(i, j - 1) = -v(i - 1, j - 1);
                    F(i - 1, j) = 0;
                }

                else if (type(i, j) == 6)
                {
                    velocity_Y(i, j) = 0;
                    //velocity_X(i, j) = 0;
                  
                   
                    F(i, j) = 0;
                    G(i, j) = 0;

                    velocity_Y(i, j - 1) = -v(i + 1, j - 1);
                    velocity_X(i - 1, j) = -u(i - 1, j + 1);
                }
                else if (type(i, j) == 7)
                {
                    velocity_Y(i, j) = 0;
                    velocity_X(i, j) = -u(i, j+1);
                    velocity_X(i - 1, j) = 0;
                    G(i, j) = 0;
                    F(i - 1, j) = 0;
                    velocity_X(i, j) = -u(i, j + 1);

                    velocity_Y(i, j - 1) = -v(i - 1, j - 1);
                }
                else if (type(i, j) == 8)
                {
                    velocity_Y(i, j) = -v(i - 1, j);
                    velocity_Y(i, j - 1) = 0;
                    //velocity_X(i, j) = 0;
                    velocity_X(i-1, j) = 0;
                    G(i, j - 1) = 0;
                    F(i - 1, j) = 0;
                    velocity_X(i, j) = -u(i, j - 1);
                }
                else if (type(i, j) == 9)
                {
                    velocity_Y(i, j) = -v(i + 1, j);
                    velocity_Y(i, j - 1) = 0;
                   // velocity_X(i, j) = 0;
                    G(i, j - 1) = 0;
                    F(i, j) = 0;
                    velocity_X(i - 1, j) = -u(i - 1, j - 1);
                }
            }
        }
    }
}