#include <iostream>
#include <cstdlib>
#include <cassert>

#include <time.h>

#include "settings.h"

// #include "pressure_solver/pressuresolver.h"
// #include "pressure_solver/sor.h"
// #include "pressure_solver/gaussseidel.h"
// #include "discretization_storage/discretization.h"
// #include "discretization_storage/donorcell.h"
// #include "discretization_storage/centraldifferences.h"
// #include "computation/computation.h"
#include "partitioning/partitioning.h"
#include "computation/computation.h"
// #include "test_and_debug/mytestfunctions.h"
//void loadFromFile(std::string filename);
#include <mpi.h>


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
Computation computation(settings);
//Partitioning partition(settings);
computation.runSimulation(settings);

  return EXIT_SUCCESS;
}