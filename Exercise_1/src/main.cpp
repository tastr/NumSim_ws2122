#include <iostream>
#include <cstdlib>
#include "settings.h"
#include "output_writer/output_writer_text.h"
#include "output_writer/output_writer_paraview.h"
#include "pressure_solver/pressuresolver.h"
#include "pressure_solver/sor.h"
#include "pressure_solver/gaussseidel.h"
#include "discretization_storage/discretization.h"
#include "discretization_storage/donorcell.h"
#include "discretization_storage/centraldifferences.h"
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
OutputWriterText myOutputWriterText(myDiscretization);
OutputWriterParaview myOutputWriterParaview(myDiscretization);

// initialize time
double current_time=0;
 //write after initialization
// myOutputWriterParaview.writeFile(current_time);
myOutputWriterText.writeFile(current_time);

// myOutputWriterParaview.writeFile(current_time);
myOutputWriterText.writeFile(current_time);

while (current_time<settings.endTime)
{
  myDiscretization->setBorderVelocity(settings.dirichletBcTop, settings.dirichletBcLeft, settings.dirichletBcRight, settings.dirichletBcBottom);
  myDiscretization->updateBoundaryFG();
  myDiscretization->updateDeltaT();
  current_time+=myDiscretization->getDeltaT();
  myDiscretization->calculation();
  myPressureSolver->calculateRHS();
  myPressureSolver->calculateP();
  myDiscretization->updateVelocity();
  
  // myOutputWriterParaview.writeFile(current_time);
  myOutputWriterText.writeFile(current_time);
}




  return EXIT_SUCCESS;
}
