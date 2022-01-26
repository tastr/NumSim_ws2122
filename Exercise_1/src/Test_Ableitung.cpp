#include <iostream>
#include <cstdlib>
#include "settings.h"
#include "output_writer/output_writer_text.h"
#include "output_writer/output_writer_paraview.h"
#include "pressure_solver/pressuresolver.h"
#include "discretization_storage/discretization.h"
#include "discretization_storage/donorcell.h"
#include "discretization_storage/centraldifferences.h"
#include <math.h>
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

CentralDifferences myCentdiff(settings);
DonorCell myDonCell(settings);


for (int j = 0; j < myCentdiff.getSize()[1]; j++)
{
  for (int i = 0; i < myCentdiff.getSize()[0]; i++)
  {
    myCentdiff.setV(i, j, sin(myCentdiff.dx() *i + myCentdiff.dy()*j));
    myCentdiff.setU(i, j, sin(myCentdiff.dx() *i + myCentdiff.dy()*j));
   
    myDonCell.setV(i, j, sin(myDonCell.dx() *i + myDonCell.dy()*j));
    myDonCell.setU(i, j, sin(myDonCell.dx() *i + myDonCell.dy()*j));
  }
  
}

     double  centdux2_sum=   0; 
     double  centduy2_sum=   0; 
     double  centdvx2_sum=   0;
     double  centdvy2_sum=   0;

     double  centdu2x_sum=   0;
     double  centdv2y_sum=   0;

     double  centduvx_sum=   0;
     double  centduvy_sum=   0;


     double  centdux_sum=   0;
     double  centdvy_sum=   0;

     double  dondux2_sum=   0 ; 
     double  donduy2_sum=   0 ; 
     double  dondvx2_sum=   0;
     double  dondvy2_sum=   0;

     double  dondu2x_sum=   0;
     double  dondv2y_sum=   0;

     double  donduvx_sum=   0;
     double  donduvy_sum=   0;

     double  dondux_sum=    0;
     double  dondvy_sum=    0;


     double  centdux2; 
     double  centduy2;
     double  centdvx2;
     double  centdvy2;

     double  centdu2x;
     double  centdv2y;

     double  centduvx;
     double  centduvy;


     double  centdux;
     double  centdvy;

     double  dondux2; 
     double  donduy2; 
     double  dondvx2;
     double  dondvy2;

     double  dondu2x;
     double  dondv2y;

     double  donduvx;
     double  donduvy;

     double  dondux;
     double  dondvy;

    double x ;
    double y ;  
    double dux2_ana;
    double duy2_ana;
    double dvx2_ana;
    double dvy2_ana;
    double du2x_ana;
    double dv2y_ana;
    double duvx_ana;
    double duvy_ana;
    double dux_ana ;
    double dvy_ana ;



for (int j = 1; j < myCentdiff.getSize()[1]-1; j++)
{
  for (int  i = 1; i < myCentdiff.getSize()[0]-1; i++)
  {
       centdux2=   myCentdiff.computeDuDx2(i,j) ; 
       centduy2=   myCentdiff.computeDuDy2(i,j) ; 
       centdvx2=   myCentdiff.computeDvDx2(i,j);
       centdvy2=   myCentdiff.computeDvDy2(i,j);

       centdu2x=   myCentdiff.computeDu2Dx(i,j);
       centdv2y=   myCentdiff.computeDv2Dy(i,j);

       centduvx=   myCentdiff.computeDuvDx(i,j);
       centduvy=   myCentdiff.computeDuvDy(i,j);


       centdux=   myDonCell.computeDuDx(i,j);
       centdvy=   myDonCell.computeDvDy(i,j);

       dondux2=   myDonCell.computeDuDx2(i,j) ; 
       donduy2=   myDonCell.computeDuDy2(i,j) ; 
       dondvx2=   myDonCell.computeDvDx2(i,j);
       dondvy2=   myDonCell.computeDvDy2(i,j);

       dondu2x=   myDonCell.computeDu2Dx(i,j);
       dondv2y=   myDonCell.computeDv2Dy(i,j);

       donduvx=   myDonCell.computeDuvDx(i,j);
       donduvy=   myDonCell.computeDuvDy(i,j);

       dondux=    myDonCell.computeDuDx(i,j);
       dondvy=    myDonCell.computeDvDy(i,j);


     x = myCentdiff.dx() *i;
     y = myCentdiff.dy()*j;  
     dux2_ana= -sin(x + y);
     duy2_ana= -sin(x + y);
     dux2_ana= -sin(x + y);
     duy2_ana= -sin(x + y);

     du2x_ana= 2*sin(x+y)*cos(x+y);
     dv2y_ana= 2*sin(x+y)*cos(x+y);
     duvx_ana= 2*sin(x+y)*cos(x+y);
     duvy_ana= 2*sin(x+y)*cos(x+y);
     
     dux_ana = cos(x+y);
     dvy_ana = cos(x+y);

       centdux2_sum=centdux2_sum + centdux2 - dux2_ana ; 
       centduy2_sum=centduy2_sum + centduy2 - duy2_ana ; 
       centdvx2_sum=centdvx2_sum + centdvx2 - dux2_ana;
       centdvy2_sum=centdvy2_sum + centdvy2 - duy2_ana;

       centdu2x_sum=centdu2x_sum + centdu2x - du2x_ana;
       centdv2y_sum=centdv2y_sum + centdv2y - dv2y_ana;

       centduvx_sum=centduvx_sum + centduvx - duvx_ana;
       centduvy_sum=centduvy_sum + centduvy - duvy_ana;

       centdux_sum=centdux_sum + centdux - dux_ana;
       centdvy_sum=centdvy_sum + centdvy - dvy_ana ;

       dondux2_sum=dondux2_sum + dondux2 - dux2_ana ; 
       donduy2_sum=donduy2_sum + donduy2 - duy2_ana  ; 
       dondvx2_sum=dondvx2_sum + dondvx2 - dux2_ana;
       dondvy2_sum=dondvy2_sum + dondvy2 - duy2_ana ;

       dondu2x_sum=dondu2x_sum + dondu2x - du2x_ana ;
       dondv2y_sum=dondv2y_sum + dondv2y - dv2y_ana;

       donduvx_sum=donduvx_sum + donduvx - duvx_ana;
       donduvy_sum=donduvy_sum + donduvy - duvy_ana;

       dondux_sum=dondux_sum + dondux - dux_ana ;
       dondvy_sum=dondvy_sum + dondvy - dvy_ana ;

     }
  }
std::cout<< "                          "  <<std::endl;
std::cout<< "Ableitungen 2. Ordnung Cent"  <<std::endl; 
std::cout<< centdux2_sum / ((myCentdiff.getSize()[1]-2)*(myCentdiff.getSize()[1]-2)) <<std::endl;
std::cout<< centduy2_sum / ((myCentdiff.getSize()[1]-2)*(myCentdiff.getSize()[1]-2)) <<std::endl;
std::cout<< centdvx2_sum / ((myCentdiff.getSize()[1]-2)*(myCentdiff.getSize()[1]-2)) <<std::endl;
std::cout<< centdvy2_sum / ((myCentdiff.getSize()[1]-2)*(myCentdiff.getSize()[1]-2)) <<std::endl;

std::cout<< "                          "  <<std::endl;
std::cout<< "Ableitungen Quadraterme Cent"  <<std::endl;
std::cout<< centdu2x_sum / ((myCentdiff.getSize()[1]-2)*(myCentdiff.getSize()[1]-2)) <<std::endl;
std::cout<< centdv2y_sum / ((myCentdiff.getSize()[1]-2)*(myCentdiff.getSize()[1]-2)) <<std::endl;

std::cout<< "                          "  <<std::endl;
std::cout<< "Ableitungen 2. uv Cent"  <<std::endl;      
std::cout<< centduvx_sum / ((myCentdiff.getSize()[1]-2)*(myCentdiff.getSize()[1]-2)) <<std::endl;
std::cout<< centduvy_sum / ((myCentdiff.getSize()[1]-2)*(myCentdiff.getSize()[1]-2)) <<std::endl;

std::cout<< "                          "  <<std::endl;
std::cout<< "Ableitungen erster ordnung Cent"  <<std::endl;
std::cout<< centdux_sum / ((myCentdiff.getSize()[1]-2)*(myCentdiff.getSize()[1]-2)) <<std::endl;
std::cout<< centdvy_sum / ((myCentdiff.getSize()[1]-2)*(myCentdiff.getSize()[1]-2)) <<std::endl;


std::cout<< "                          "  <<std::endl;
std::cout<< "Ableitungen 2. Ordnung Don"  <<std::endl; 
std::cout<< dondux2_sum / ((myCentdiff.getSize()[1]-2)*(myCentdiff.getSize()[1]-2)) <<std::endl;
std::cout<< donduy2_sum / ((myCentdiff.getSize()[1]-2)*(myCentdiff.getSize()[1]-2)) <<std::endl;
std::cout<< dondvx2_sum / ((myCentdiff.getSize()[1]-2)*(myCentdiff.getSize()[1]-2)) <<std::endl;
std::cout<< dondvy2_sum / ((myCentdiff.getSize()[1]-2)*(myCentdiff.getSize()[1]-2)) <<std::endl;

std::cout<< "                          "  <<std::endl;
std::cout<< "Ableitungen Quadraterme Don"  <<std::endl;
std::cout<< dondu2x_sum / ((myCentdiff.getSize()[1]-2)*(myCentdiff.getSize()[1]-2)) <<std::endl;
std::cout<< dondv2y_sum / ((myCentdiff.getSize()[1]-2)*(myCentdiff.getSize()[1]-2)) <<std::endl;

std::cout<< "                          "  <<std::endl;
std::cout<< "Ableitungen uv Don"  <<std::endl; 
std::cout<< donduvx_sum / ((myCentdiff.getSize()[1]-2)*(myCentdiff.getSize()[1]-2)) <<std::endl;
std::cout<< donduvy_sum / ((myCentdiff.getSize()[1]-2)*(myCentdiff.getSize()[1]-2)) <<std::endl;

std::cout<< "                          "  <<std::endl;
std::cout<< "Ableitungen 1. Ordnung Don"  <<std::endl;      
std::cout<< dondux_sum / ((myCentdiff.getSize()[1]-2)*(myCentdiff.getSize()[1]-2)) <<std::endl;
std::cout<< dondvy_sum / ((myCentdiff.getSize()[1]-2)*(myCentdiff.getSize()[1]-2)) <<std::endl;
std::cout<< "                          "  <<std::endl;



  return EXIT_SUCCESS;
}
