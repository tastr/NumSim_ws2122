#include <iostream>
#include <cstdlib>
#include <array>
#include <fstream>
#include <iomanip>
#include "array2d.h"
#include "fieldvariable.h"
#include "staggeredgrid.h"
#include "settings.h"
#include "discretization.h"
#include "sor.h"
#include "pressuresolver.h"
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
for (size_t i = 0; i < 6; i++)
{
mydiscretization.pressure(i,i)=i;
mydiscretization.calculation();
mydiscretization.pressure(i,2)=i;
mydiscretization.pressure(i,3)=i;
mydiscretization.pressure(2,i)=i;
mydiscretization.pressure(4,i)=i;
mydiscretization.pressure(3,i)=i;

mydiscretization.F(i,i)=i;
mydiscretization.calculation();
mydiscretization.F(i,2)=i;
mydiscretization.F(i,3)=i;
mydiscretization.F(2,i)=i;
mydiscretization.F(4,i)=i;
mydiscretization.F(3,i)=i;

}


mydiscretization.print("F");

//PressureSolver mypressuresolver(mydiscretization);
//mypressuresolver.setPressureBoundaries();

mydiscretization.pressure(1,1)=1;
mydiscretization.pressure(2,1)=4;
mydiscretization.pressure(4,4)=2;
mydiscretization.print("pressure");
std::cout <<mydiscretization.pressure(5,5) << std::endl;
SOR mysor(mydiscretization);
mysor.Iterationsverfahren();



return EXIT_SUCCESS;
}
