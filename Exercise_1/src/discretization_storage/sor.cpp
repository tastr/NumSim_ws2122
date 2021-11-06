#include "array2d.h"
#include "pressuresolver.h"
#include "discretization.h"
#include "sor.h"
#include <cassert>

SOR::SOR(Discretization discretization_)
:PressureSolver(discretization_),discretization(discretization_)
{
    float tol = 0.000001;

}

SOR::~SOR()
{
}


//FieldVariable SOR::Iterationsverfahren()
void SOR::Iterationsverfahren()
    {
    FieldVariable p   = discretization.p();
    float deltax_quad = discretization.dx() * discretization.dx();
    float deltay_quad = discretization.dy() * discretization.dy();
    float vorfaktor= deltax_quad * deltay_quad/ (2 * (deltay_quad + deltax_quad ));
    int i_max = discretization.getSize()[0] , j_max = discretization.getSize()[1];
    int safe=0;
    do
    { FieldVariable p  = discretization.p();
       for (int j = 1; j < j_max-1 ; j++)
        {   for (int i = 1; i < i_max-1 ; i++)
            {   
            p(i,j) = vorfaktor*( (p(i-1,j) + p(i+1,j))/deltax_quad + (p(i,j-1) + p(i,j+1)) / deltay_quad   - discretization.rhs(i,j)) ;
            //residuum = residuum + abs_(p(i,j)-discretization.RHS(i,j));   
            //residuum =  residuum +abs_( (p(i-1,j)  - 2 * p(i,j) + p(i+1,j))/deltax_quad + (p(i,j-1) - 2 * p(i,j)  + p(i,j+1)) / deltay_quad   - discretization.RHS(i,j));
            }            
        } 
        setPressureBoundaries();
        discretization.setP(p) ;     
        safe++;
        std::cout<< "Residuum " << residuum() << " Safe "<< safe <<std::endl;
    }while(residuum() > 0.00001 && safe<1000);
    //p.print();
    }


