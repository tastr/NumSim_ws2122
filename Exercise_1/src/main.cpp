#include <iostream>
#include <cstdlib>
#include "settings.h"
#include "output_writer/write_paraview_output.h"

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

  // write 5 output files
  for (int i = 0; i < 5; i++)
  {
    writeParaviewOutput(i);
  }
  return EXIT_SUCCESS;
}
