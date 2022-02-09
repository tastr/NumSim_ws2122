#include "pressuresolver.h"
#include <cassert>

PressureSolver::PressureSolver(Discretization &discretization)
    : discretization_(discretization)
{
  settings_ = discretization_.getSettings();
}

PressureSolver::~PressureSolver()
{
}

double PressureSolver::abs_(double number)
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

void PressureSolver::setPressureBoundaries()
{
  int i_max = discretization_.pIEnd(), j_max = discretization_.pJEnd();
  int i_begin = discretization_.pIBegin(), j_begin = discretization_.pJBegin();

  /*
  for (int j = j_begin; j < j_max; j++)
  {
    discretization_.setP(i_begin, j, discretization_.p(i_begin + 1, j));
    discretization_.setP(i_max - 1, j, discretization_.p(i_max - 2, j));
  }

  for (int i = i_begin; i < i_max; i++)
  {
    discretization_.setP(i, j_begin, discretization_.p(i, j_begin + 1));
    discretization_.setP(i, j_max - 1, discretization_.p(i, j_max - 2));
  }
  */

  if (settings_.outflowLeft == true)
    {
      for (int j = j_begin; j < j_max; j++)
      {
        discretization_.setP(i_begin, j, 0);
      }
    } 
    else
    {
      for (int j = j_begin; j < j_max; j++)
      {
        discretization_.setP(i_begin, j, discretization_.p(i_begin + 1, j));
      }
    }

    if (settings_.outflowRight == true)
    {
      for (int j = j_begin; j < j_max; j++)
      {
        discretization_.setP(i_max - 1, j, 0);
      }
    } 
    else
    {
      for (int j = j_begin; j < j_max; j++)
      {
        discretization_.setP(i_max - 1, j, discretization_.p(i_max - 2, j));
      }
    }

    if (settings_.outflowBottom == true)
    {
      for (int i = i_begin; i < i_max; i++)
      {
        discretization_.setP(i, j_begin, 0);
      }
    } 
    else
    {
      for (int i = i_begin; i < i_max; i++)
      {
        discretization_.setP(i, j_begin, discretization_.p(i, j_begin + 1));
      }
    }

    if (settings_.outflowTop == true)
    {
      for (int i = i_begin; i < i_max; i++)
      {
        discretization_.setP(i, j_max - 1, 0);
      }
    } 
    else
    {
      for (int i = i_begin; i < i_max; i++)
      {
        discretization_.setP(i, j_max - 1, discretization_.p(i, j_max - 2));
      }
    }

}

void PressureSolver::calculateRHS()
{
  double Fij = 0;
  double Fim1j = 0;
  double Gij = 0;
  double Gijm1;
  double current_rhs = 0;
  double sumRHS = 0;
  for (int j = discretization_.pJBegin() + 1; j < discretization_.pJEnd() - 1; j++)
  {
    for (int i = discretization_.pIBegin() + 1; i < discretization_.pIEnd() - 1; i++)
    {
      if (discretization_.getTyp(i, j) == 0)
      {
        Fij = discretization_.f(i, j);
        Fim1j = discretization_.f(i - 1, j);
        Gij = discretization_.g(i, j);
        Gijm1 = discretization_.g(i, j - 1);

        current_rhs = ((Fij - Fim1j) / discretization_.dx() + (Gij - Gijm1) / discretization_.dy()) / discretization_.getDeltaT();
        current_rhs = ( 1- discretization_.getURrhs()) * discretization_.rhs(i, j) + discretization_.getURrhs() * current_rhs;
        discretization_.setRHS(i, j, current_rhs);
        sumRHS += current_rhs;
      }
    }
  }
  std::cout << "Summe RHS: " << sumRHS << "\t";
}

void PressureSolver::calculateP()
{
  // is virtual class only should never be called
  assert(false);
}

double PressureSolver::residuum()
{
  double res = 0;
  double pxx = 0;
  double pyy = 0;
  double temp= 0;
  for (int j = discretization_.pJBegin() + 1; j < discretization_.pJEnd() - 1; j++)
  {
    for (int i = discretization_.pIBegin() + 1; i < discretization_.pIEnd() - 1; i++)
    {
      if (discretization_.getTyp(i, j) == 0)
      {
      pxx = (discretization_.p(i - 1, j) - 2 * discretization_.p(i, j) + discretization_.p(i + 1, j)) / (discretization_.dx() * discretization_.dx());
      pyy = (discretization_.p(i, j - 1) - 2 * discretization_.p(i, j) + discretization_.p(i, j + 1)) / (discretization_.dy() * discretization_.dy());
      temp = pxx + pyy - discretization_.rhs(i, j);
      res = res + (temp * temp);
      //printf("%2.5f ",1000*temp);
      }else
      {
      //printf("9.99999");
      }
    }
    //printf("\n");
  }
   //printf("\n\n\n\n");
  return res;
}
