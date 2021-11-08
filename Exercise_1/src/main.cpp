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
myOutputWriterParaview.writeFile(current_time);
myOutputWriterText.writeFile(current_time);
myDiscretization->setBorderVelocity(settings.dirichletBcTop, settings.dirichletBcLeft, settings.dirichletBcRight, settings.dirichletBcBottom);
myOutputWriterParaview.writeFile(current_time);
myOutputWriterText.writeFile(current_time);

// while (current_time<settings.endTime)
// {
  myDiscretization->updateDeltaT();
  current_time+=myDiscretization->getDeltaT();
  myDiscretization->calculation();
  myPressureSolver->calculateRHS();
  myPressureSolver->calculateP();
  
  
  myOutputWriterParaview.writeFile(current_time);
  myOutputWriterText.writeFile(current_time);
// }





// Declare an object of type donorcell.
//   DonorCell mydonorcell;
   // Declare two pointers, one of type DonorCell * and the other
   //  of type Dis *, and initialize them to point to mydonorcell.
  // DonorCell *p_mydonorcell = &mydonorcell;
   //Discretization    *p_discretization = &mydonorcell;

   // Call the functions.
   //pdiscretization->function();           // Call virtual function.
   //pdonorcell->function(); // Call nonvirtual function.


// std::shared_ptr<CentralDifferences> centralDifferences1 = std::make_shared<CentralDifferences>(settings);
// double value1 = centralDifferences1->computeDuvDx(1,1);

// std::shared_ptr<Discretization> centralDifferences2 = std::make_shared<CentralDifferences>(settings);
// double value2 = centralDifferences2->computeDuvDx(1,1);





  return EXIT_SUCCESS;
}
