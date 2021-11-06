#include <iostream>
#include <cstdlib>
#include "settings.h"
#include "output_writer/output_writer_text.h"
#include "output_writer/output_writer_paraview.h"
#include "pressure_solver/pressuresolver.h"
#include "discretization_storage/discretization.h"

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
Discretization myDiscretization(settings);
PressureSolver myPressureSolver(myDiscretization); //TODO use reference instead
OutputWriterText myOutputWriterText(std::make_shared<Discretization>(myDiscretization));
OutputWriterParaview myOutputWriterParaview(std::make_shared<Discretization>(myDiscretization));

// initialize time
double current_time=0;
// write after initialization
myOutputWriterParaview.writeFile(current_time);
myOutputWriterText.writeFile(current_time);
myDiscretization.setBorderVelocity(settings.dirichletBcTop, settings.dirichletBcLeft, settings.dirichletBcRight, settings.dirichletBcBottom);
myOutputWriterParaview.writeFile(current_time);
myOutputWriterText.writeFile(current_time);

myDiscretization.updateDeltaT();
current_time+=myDiscretization.getDeltaT();
myDiscretization.calculation();

myOutputWriterParaview.writeFile(current_time);
myOutputWriterText.writeFile(current_time);


  return EXIT_SUCCESS;
}
