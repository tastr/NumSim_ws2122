#include "sor.h"
#include <cassert>

SOR::SOR(Discretization &discretization_)
    : PressureSolver(discretization_)
{
}

SOR::~SOR()
{
}

void SOR::calculateP()
{
  double deltax_quad = discretization_.dx() * discretization_.dx();
  double deltay_quad = discretization_.dy() * discretization_.dy();
  double vorfaktor = deltax_quad * deltay_quad / (2 * (deltay_quad + deltax_quad));
  int i_max = discretization_.getSize()[0], j_max = discretization_.getSize()[1];
  int safe = 0;

  // auxiliarys for code readability
  double x_term;
  double y_term;
  double omega = discretization_.getOmega();

  double epsilonquad = discretization_.getepsilon() * discretization_.getepsilon();
  double resterm;
  int Nnumber = (discretization_.nCells()[0] * discretization_.nCells()[1]);

  do
  {

    for (int j = 1; j < j_max - 1; j++)
    {
      for (int i = 1; i < i_max - 1; i++)
      {
        if (discretization_.getTyp(i, j) == 0)
        {
          x_term = (discretization_.p(i - 1, j) + discretization_.p(i + 1, j)) / deltax_quad;
          y_term = (discretization_.p(i, j - 1) + discretization_.p(i, j + 1)) / deltay_quad;
          // value of pij gets overwritten with the new approximation
          discretization_.setP(i, j, (1 - omega) * discretization_.p(i, j) + omega * vorfaktor * (x_term + y_term - discretization_.rhs(i, j)));

          // Gaussseideltermijend= v*((discretization_.p(i-1,jend)+discretization_.p(i+1,jend))/dx2+ discretization_.p(i,jend-1)/dy2-discretization_.rhs(i,jend))/div_y;
          // discretization_.setP(i,jend, (1-omega) * discretization_.p(i,jend)+ omega* Gaussseideltermijend);
        }
      }
    }
    safe++;
    discretization_.updatedPressureBC();
    discretization_.setObstaclePressure();
    resterm = (residuum() * residuum()) / Nnumber;

  } while (resterm > epsilonquad && safe < 20000);
  // std::cout<< "Residuum " << residuum() << " Safe "<< safe <<std::endl;
  // setPressureBoundaries();
  // discretization_.updatedPressureBC();
}
