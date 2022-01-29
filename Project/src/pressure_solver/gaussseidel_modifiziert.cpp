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

  double dx2 = deltax_quad;
  double dy2 = deltay_quad;
  int i_max = discretization_.getSize()[0], j_max = discretization_.getSize()[1];

  int iend = i_max - 2;
  int jend = j_max - 2;
  double div_ecke = 1.0 - vorfaktor * (1.0 / dx2 + 1.0 / dy2);
  double div_x = 1.0 - vorfaktor * (1.0 / dx2);
  double div_y = 1.0 - vorfaktor * (1.0 / dy2);
  double v = vorfaktor;

  do
  { // FieldVariable p=discretization_.p();

    discretization_.setP(1, 1, v * (discretization_.p(2, 1) / dx2 + discretization_.p(1, 2) / dy2 - discretization_.rhs(1, 1)) / div_ecke);
    discretization_.setP(1, jend, v * (discretization_.p(2, jend) / dx2 + discretization_.p(1, jend - 1) / dy2 - discretization_.rhs(1, jend)) / div_ecke);
    discretization_.setP(iend, 1, v * (discretization_.p(iend - 1, 1) / dx2 + discretization_.p(iend, 2) / dy2 - discretization_.rhs(iend, 1)) / div_ecke);
    discretization_.setP(iend, jend, v * (discretization_.p(iend - 1, jend) / dx2 + discretization_.p(iend, jend - 1) / dy2 - discretization_.rhs(iend, jend)) / div_ecke);

    for (int j = discretization_.pJBegin() + 2; j < discretization_.pJEnd() - 2; j++)

    {
      discretization_.setP(1, j, v * (discretization_.p(2, j) / dx2 + (discretization_.p(1, j - 1) + discretization_.p(1, j + 1)) / dy2 - discretization_.rhs(1, j)) / div_x);
      discretization_.setP(j, 1, v * ((discretization_.p(j - 1, 1) + discretization_.p(j + 1, 1)) / dx2 + discretization_.p(j, 2) / dy2 - discretization_.rhs(j, 1)) / div_y);

      for (int i = discretization_.pIBegin() + 2; i < discretization_.pIEnd() - 2; i++)
      {
        x_term = (discretization_.p(i - 1, j) + discretization_.p(i + 1, j)) / deltax_quad;
        y_term = (discretization_.p(i, j - 1) + discretization_.p(i, j + 1)) / deltay_quad;
        // value of pij gets overwritten with the new approximation

        discretization_.setP(i, j, vorfaktor * (x_term + y_term - discretization_.rhs(i, j)));

        // x_term= (p(i-1,j) + p(i+1,j))/deltax_quad;
        // y_term= (p(i,j-1) + p(i,j+1)) / deltay_quad;

        // discretization_.p(i,j)=vorfaktor*( x_term + y_term   - discretization_.rhs(i,j));
      }

      discretization_.setP(iend, j, v * (discretization_.p(iend - 1, j) / dx2 + (discretization_.p(iend, j - 1) + discretization_.p(iend, j + 1)) / dy2 - discretization_.rhs(iend, j)) / div_x);
      discretization_.setP(j, jend, v * ((discretization_.p(j - 1, jend) + discretization_.p(j + 1, jend)) / dx2 + discretization_.p(j, jend - 1) / dy2 - discretization_.rhs(j, jend)) / div_y);
    }
    safe++;

    // discretization_.setP(p);
    // Ohne den konvergiert er doppelt so schnell
    discretization_.updatedPressureBC();

  } while (residuum() > discretization_.getepsilon() && safe < 20000);
  std::cout << "Residuum " << residuum() << " Safe " << safe << std::endl;
  discretization_.updatedPressureBC();
}
