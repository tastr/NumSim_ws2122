#include "gaussseidel.h"
#include <cassert>

GaussSeidel::GaussSeidel(Discretization &discretization_)
    : PressureSolver(discretization_)
{
}

GaussSeidel::~GaussSeidel()
{
}

void GaussSeidel::calculateP()
{ // FieldVariable p=discretization_.p();

  double deltax_quad = discretization_.dx() * discretization_.dx();
  double deltay_quad = discretization_.dy() * discretization_.dy();
  double vorfaktor = deltax_quad * deltay_quad / (2 * (deltay_quad + deltax_quad));
  // FieldVariable p=discretization_.p();
  int safe = 0;

  // auxilliarys for code readability
  double x_term;
  double y_term;
  double omega = discretization_.getOmega();

  double epsilonquad = discretization_.getepsilon() * discretization_.getepsilon();
  double resterm;
  int Nnumber = (discretization_.nCells()[0] * discretization_.nCells()[1]);

  do
  {

    for (int j = discretization_.pJBegin() + 1; j < discretization_.pJEnd() - 1; j++)

    {

      for (int i = discretization_.pIBegin() + 1; i < discretization_.pIEnd() - 1; i++)
      {
        if (discretization_.getTyp(i, j) == 0)
        {
          x_term = (discretization_.p(i - 1, j) + discretization_.p(i + 1, j)) / deltax_quad;
          y_term = (discretization_.p(i, j - 1) + discretization_.p(i, j + 1)) / deltay_quad;

          discretization_.setP(i, j, vorfaktor * (x_term + y_term - discretization_.rhs(i, j)));
        }
      }
    }
    safe++;

    discretization_.updatedPressureBC();
    discretization_.setObstaclePressure();
    resterm = residuum() / Nnumber;

  } while (resterm > epsilonquad && safe < 20000);
  std::cout << "Residuum " << residuum() << " Safe " << safe << std::endl;
}
