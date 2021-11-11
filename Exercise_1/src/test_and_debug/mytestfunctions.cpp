#include "mytestfunctions.h"
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <memory>

bool MyTestFunctions::testFunction1(Settings settings)
{
  //create objects of classes
std::shared_ptr<Discretization> myTestDiscretization;
if (settings.useDonorCell)
{
  myTestDiscretization=std::make_shared<DonorCell>(settings);
} else
{
  myTestDiscretization=std::make_shared<CentralDifferences>(settings);
}
std::shared_ptr<PressureSolver> myTestPressureSolver;
if (settings.pressureSolver == "SOR")
{
  myTestPressureSolver= std::make_shared<SOR>(*myTestDiscretization); //reference is used
} else
{
  myTestPressureSolver= std::make_shared<GaussSeidel>(*myTestDiscretization); //reference is used
}
OutputWriterText myOutputWriterText(myTestDiscretization);
OutputWriterParaview myOutputWriterParaview(myTestDiscretization);

double p_value=99;
double current_time=0;

myTestDiscretization->setBorderVelocity(settings.dirichletBcTop, settings.dirichletBcLeft, settings.dirichletBcRight, settings.dirichletBcBottom);
myTestDiscretization->updateBoundaryFG();

myTestDiscretization->updateDeltaT();
current_time+=myTestDiscretization->getDeltaT();

for (int j = 1; j < settings.nCells[1]+1; j++)
{
  for (int i = 1; i < settings.nCells[0]+1; i++)
  {
    myTestDiscretization->setP(i,j,p_value);
  }
  
}

myTestPressureSolver->setPressureBoundaries();
myOutputWriterText.writeFile(current_time);

myTestPressureSolver->calculateP();
myOutputWriterText.writeFile(current_time);


// Test2
double rhs_value= 100;
p_value=0;
for (int j = 1; j < settings.nCells[1]+1; j++)
{
  for (int i = 1; i < settings.nCells[0]+1; i++)
  {
    myTestDiscretization->setRHS(i,j,rhs_value);
    myTestDiscretization->setP(i,j,p_value);
  }
  
}

myTestPressureSolver->setPressureBoundaries();
myOutputWriterText.writeFile(current_time);

myTestPressureSolver->calculateP();
myOutputWriterText.writeFile(current_time);

  return true;
}

bool MyTestFunctions::testFunction2(Settings settings)
{
    //create objects of classes
  std::shared_ptr<Discretization> myTestDiscretization;
  if (settings.useDonorCell)
  {
    myTestDiscretization=std::make_shared<DonorCell>(settings);
  } else
  {
    myTestDiscretization=std::make_shared<CentralDifferences>(settings);
  }
  std::shared_ptr<PressureSolver> myTestPressureSolver;
  if (settings.pressureSolver == "SOR")
  {
    myTestPressureSolver= std::make_shared<SOR>(*myTestDiscretization); //reference is used
  } else
  {
    myTestPressureSolver= std::make_shared<GaussSeidel>(*myTestDiscretization); //reference is used
  }
  OutputWriterText myOutputWriterText(myTestDiscretization);
  OutputWriterParaview myOutputWriterParaview(myTestDiscretization);

  double p_value=99;
  double current_time=0;

  myTestDiscretization->setBorderVelocity(settings.dirichletBcTop, settings.dirichletBcLeft, settings.dirichletBcRight, settings.dirichletBcBottom);
  myTestDiscretization->updateBoundaryFG();

  myTestDiscretization->updateDeltaT();
  current_time+=myTestDiscretization->getDeltaT();


  // Test2
  double rhs_value= 1;
  p_value=0;
  for (int j = 1; j < settings.nCells[1]+1; j++)
  {
    for (int i = 1; i < settings.nCells[0]+1; i++)
    {
      // p_value=i*i+j*j;
      myTestDiscretization->setRHS(i,j,(rhs_value+i));
      myTestDiscretization->setP(i,j,p_value);
    }
    
  }
  // myTestDiscretization->setRHS(1,1,6);
  // myTestDiscretization->setRHS(1,2,0);
  // myTestDiscretization->setRHS(2,1,0);
  // myTestDiscretization->setRHS(2,2,-6);
  // myTestDiscretization->setP(1,1,0);

  myTestPressureSolver->setPressureBoundaries();
  myOutputWriterText.writeFile(current_time);

  myTestPressureSolver->calculateP();
  myOutputWriterText.writeFile(current_time);
  

    return true;
}