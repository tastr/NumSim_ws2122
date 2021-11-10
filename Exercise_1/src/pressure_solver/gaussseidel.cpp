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
    int safe=0;
    do
    { 
       for (int j = 1; j < j_max-1; j++)
        {
       for (int i = 1; i < i_max-1; i++)
            {   
            discretization_.setP(i,j, vorfaktor*( (discretization_.p(i-1,j) + discretization_.p(i+1,j))/deltax_quad + (discretization_.p(i,j-1) + discretization_.p(i,j+1)) / deltay_quad   - discretization_.rhs(i,j))) ;
            }  
                 
        }        
         safe++;
        
    }while(residuum() > discretization_.getepsilon()  && safe<2000);
    std::cout<< "Residuum " << residuum() << " Safe "<< safe <<std::endl;
    discretization_.updatedPressureBC();
    }

    
