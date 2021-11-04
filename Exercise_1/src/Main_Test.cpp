#include "discretization_storage/array2d.h"
#include "discretization_storage/fieldvariable.h"
#include "discretization_storage/staggeredgrid.h"
#include "settings.h"
#include "discretization_storage/discretization.h"
#include "output_writer/output_writer.h"

#include <iostream>
#include <cstdlib>
#include <array>
#include <fstream>
#include <iomanip>

//g++ fieldvariable.cpp staggeredgrid.cpp  array2d.cpp Main_Test.cpp settings.cpp -o testA

//#include "centraldifferences.h"

int main(int argc, char *argv[]){
	

std::array<int,2> size{2,2};
Array2D array2D(size);        // object of class Array2D with size 2x2

array2D(0,1) = 1.0;
double value = array2D(0,1);
std::cout << value << std::endl;

std::cout<<"Diese Funktion printed die Fieldvariable" <<std::endl;
size={5,5};
FieldVariable Test(size);
Test.print();
Test(0,2)=5;
std::cout<< "" <<std::endl;
Test.print();
std::cout<< "Maximum " << Test.max() <<std::endl;

std::array<double,2> top={1,2};
std::array<double,2> left={3,4};
std::array<double,2> right={5,6};
std::array<double,2> bottom={7,8};

StaggeredGrid testgrid(size);
testgrid.setBorderVelocity(top,left,right,bottom);

testgrid.print("velocity_x");
std::cout<< " " <<std::endl;
testgrid.print("velocity_y");



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

Discretization mydiscretization(settings);
mydiscretization.print("velocity_x");
//mydiscretization.print("velocity_x");

mydiscretization.calculation();
mydiscretization.print("F");

return EXIT_SUCCESS;
}
