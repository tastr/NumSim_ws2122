#include <iostream>
#include <cstdlib>
#include <cassert>

#include <time.h>

#include "settings.h"
#include "output_writer/output_writer_text.h"
#include "output_writer/output_writer_paraview.h"
#include "pressure_solver/pressuresolver.h"
#include "pressure_solver/sor.h"
#include "pressure_solver/gaussseidel.h"
#include "discretization_storage/discretization.h"
#include "discretization_storage/donorcell.h"
#include "discretization_storage/centraldifferences.h"
#include "test_and_debug/mytestfunctions.h"
//void loadFromFile(std::string filename);



int main(int argc, char *argv[])
{
  // if the number of given command line arguments is only 1 (= the program name), print out usage information and exit
  if (argc == 1)
  {
    std::cout << "usage: " << argv[0] << " <filename>" << std::endl;

    return EXIT_FAILURE;
  }

  // read in the first argument
  std::string filename = argv[1];

  // print message
  std::cout << "Filename: \"" << filename << "\"" << std::endl;
Settings settings;
// load settings from file
settings.loadFromFile(filename);
// display all settings on console
settings.printSettings();

// testing 
// MyTestFunctions myTest;
// assert(myTest.testFunction2(settings));
double time1, tstart;
tstart=clock();

//create objects of classes
std::shared_ptr<Discretization> myDiscretization;
if (settings.useDonorCell)
{
  myDiscretization=std::make_shared<DonorCell>(settings);
} else
{
  myDiscretization=std::make_shared<CentralDifferences>(settings);
}
std::shared_ptr<PressureSolver> myPressureSolver;
if (settings.pressureSolver == "SOR")
{
  myPressureSolver= std::make_shared<SOR>(*myDiscretization); //reference is used
} else
{
  myPressureSolver= std::make_shared<GaussSeidel>(*myDiscretization); //reference is used
}

// std::shared_ptr<Discretization> pointer_to_myDiscretization (& myDiscretization); //vermutlich gibt es da einen besseren Weg, aber den habe ich nicht gefunden...
// OutputWriterText myOutputWriterText(myDiscretization);
OutputWriterParaview myOutputWriterParaview(myDiscretization);
int Laufzaehler=0;
// initialize time
double current_time=0;
 //write after initialization

myDiscretization->setBorderVelocity(settings.dirichletBcTop, settings.dirichletBcLeft, settings.dirichletBcRight, settings.dirichletBcBottom);
myDiscretization->updateBoundaryFG();
// myOutputWriterParaview.writeFile(current_time);
// myOutputWriterText.writeFile(current_time);

while (current_time<settings.endTime)
{
  myDiscretization->updateDeltaT();
  current_time+=myDiscretization->getDeltaT();
  myDiscretization->calculation();
  myPressureSolver->calculateRHS();
  myPressureSolver->calculateP();
  myDiscretization->updateVelocity();
  myDiscretization->setBorderVelocity(settings.dirichletBcTop, settings.dirichletBcLeft, settings.dirichletBcRight, settings.dirichletBcBottom);
  myDiscretization->updateBoundaryFG();
  
  myOutputWriterParaview.writeFile(current_time);
  // myOutputWriterText.writeFile(current_time);
Laufzaehler=Laufzaehler+1;
}
std::cout<< "Noetige Iterationen " << Laufzaehler <<std::endl;
time1=clock()-tstart;
time1=time1/CLOCKS_PER_SEC;

std::cout<< "Laufzeit in s " << time1  <<std::endl;

  return EXIT_SUCCESS;
}


