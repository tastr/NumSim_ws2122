#include "gaussseidel.h"
#include <cassert>

GaussSeidel::GaussSeidel(Discretization discretization_)
:PressureSolver(discretization_)
{
    double tol = 0.000001;
    
}

GaussSeidel::~GaussSeidel()
{
}


//FieldVariable SOR::Iterationsverfahren()
void GaussSeidel::Iterationsverfahren()
    {
    FieldVariable p   = discretization_.p();
    double deltax_quad = discretization_.dx() * discretization_.dx();
    double deltay_quad = discretization_.dy() * discretization_.dy();
    double vorfaktor= deltax_quad * deltay_quad/ (2 * (deltay_quad + deltax_quad ));
    int i_max = discretization_.getSize()[0] , j_max = discretization_.getSize()[1];
    int safe=0;
    do
    { FieldVariable p  = discretization_.p();
       for (int j = 1; j < j_max-1 ; j++)
        {   for (int i = 1; i < i_max-1 ; i++)
            {   
            p(i,j) = vorfaktor*( (p(i-1,j) + p(i+1,j))/deltax_quad + (p(i,j-1) + p(i,j+1)) / deltay_quad   - discretization_.rhs(i,j)) ;
            
            //residuum = residuum + abs_(p(i,j)-discretization_.RHS(i,j));   
            //residuum =  residuum +abs_( (p(i-1,j)  - 2 * p(i,j) + p(i+1,j))/deltax_quad + (p(i,j-1) - 2 * p(i,j)  + p(i,j+1)) / deltay_quad   - discretization.RHS(i,j));
            }            
        } 
        setPressureBoundaries();
        discretization_.setP(p) ;     
        safe++;
        std::cout<< "Residuum " << residuum() << " Safe "<< safe <<std::endl;
    }while(residuum() > 0.00001 && safe<1000);
    //p.print();
    }


