#include "gaussseidel.h"
#include <cassert>

GaussSeidel::GaussSeidel(Discretization& discretization_)
:PressureSolver(discretization_)
{    
}

GaussSeidel::~GaussSeidel()
{
}



void GaussSeidel::calculateP()
    {
    double deltax_quad = discretization_.dx() * discretization_.dx();
    double deltay_quad = discretization_.dy() * discretization_.dy();
    double vorfaktor= deltax_quad * deltay_quad/ (2 * (deltay_quad + deltax_quad ));
    int i_max = discretization_.getSize()[0] , j_max = discretization_.getSize()[1];
    // FieldVariable p=discretization_.p();
    int safe=0;

    //auxilliarys for code readability
    double x_term;
    double y_term;
    double omega=discretization_.getOmega();
    do
    { 
       for (int j = 1; j < j_max-1; j++)
        {
       for (int i = 1; i < i_max-1; i++)
            {   
            x_term= (discretization_.p(i-1,j) + discretization_.p(i+1,j))/deltax_quad;
            y_term= (discretization_.p(i,j-1) + discretization_.p(i,j+1)) / deltay_quad;  
              // value of pij gets overwritten with the new approximation
          
            discretization_.setP(i,j, vorfaktor*( x_term + y_term   - discretization_.rhs(i,j))) ;
            // p(i,j)=vorfaktor*( x_term + y_term   - discretization_.rhs(i,j));
            }  
                 
        }        
         safe++;
        //  discretization_.setP(p);
        discretization_.updatedPressureBC();
        
        
    }while(residuum() > discretization_.getepsilon()  && safe<2000);
    std::cout<< "Residuum " << residuum() << " Safe "<< safe <<std::endl;
    
    
    }

    
