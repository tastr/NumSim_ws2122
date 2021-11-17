#include "sor.h"
#include <cassert>

SOR::SOR(Discretization& discretization_)
:PressureSolver(discretization_)
{
}

SOR::~SOR()
{
}



void SOR::calculateP()
    {
    double deltax_quad = discretization_.dx() * discretization_.dx();
    double deltay_quad = discretization_.dy() * discretization_.dy();
    double vorfaktor= deltax_quad * deltay_quad/ (2 * (deltay_quad + deltax_quad ));
    int i_max = discretization_.getSize()[0] , j_max = discretization_.getSize()[1];
    int safe=0;

    
    //auxiliarys for code readability    
    double x_term;
    double y_term;
    double omega=discretization_.getOmega();
    

    double dx2=deltax_quad;
    double dy2=deltay_quad;
    
    int iend=i_max-2;
    int jend=j_max-2;
    double div_ecke=1.0-vorfaktor*(1.0/dx2+1.0/dy2);
    double div_x   =1.0-vorfaktor*(1.0/dx2);
    double div_y   =1.0-vorfaktor*(1.0/dy2);
    double v=vorfaktor;

     double  Gaussseidelterm11;                  
     double  Gaussseidelterm1jend;       
     double  Gaussseideltermiend1;  
     double  Gaussseideltermiendjend; 

    double Gaussseidelterm1j;
    double Gaussseideltermj1;
    double Gaussseideltermiendj ;
    double Gaussseideltermijend;

    

    do
    {
       Gaussseidelterm11=      v*(discretization_.p(2,1)/dx2+discretization_.p(1,2)/dy2 -discretization_.rhs(1,1))/div_ecke; 
       Gaussseidelterm1jend=   v*(discretization_.p(2,jend)/dx2+discretization_.p(1,jend-1)/dy2 -discretization_.rhs(1,jend))/div_ecke;
       Gaussseideltermiend1=   v*(discretization_.p(iend-1,1)/dx2+discretization_.p(iend,2)/dy2        -discretization_.rhs(iend,1))/div_ecke;
       Gaussseideltermiendjend=v*(discretization_.p(iend-1,jend)/dx2+discretization_.p(iend,jend-1)/dy2-discretization_.rhs(iend,jend))/div_ecke;

       
       discretization_.setP(1,1,       (1-omega) * discretization_.p(1,1)       + omega*  Gaussseidelterm11);
       discretization_.setP(1,jend,    (1-omega) * discretization_.p(1,jend)    + omega*  Gaussseidelterm1jend); 
       discretization_.setP(iend,1,    (1-omega) * discretization_.p(iend,1)    + omega*  Gaussseideltermiend1);
       discretization_.setP(iend,jend, (1-omega) * discretization_.p(iend,jend) + omega*  Gaussseideltermiendjend); 

       for (int i = 2; i < i_max-2; i++)
        {
            Gaussseideltermj1= v*((discretization_.p(i-1,1)+discretization_.p(i+1,1))/dx2+ discretization_.p(i,2)/ dy2-discretization_.rhs(i,1))/div_y;      
            discretization_.setP(i,1, (1-omega) * discretization_.p(i,1)+ omega*  Gaussseideltermj1);
            
            Gaussseideltermijend= v*((discretization_.p(i-1,jend)+discretization_.p(i+1,jend))/dx2+ discretization_.p(i,jend-1)/dy2-discretization_.rhs(i,jend))/div_y;
            discretization_.setP(i,jend, (1-omega) * discretization_.p(i,jend)+ omega* Gaussseideltermijend);
           
        }

           
       for (int j = 2; j < j_max-2 ; j++)
        {   Gaussseidelterm1j= v*(discretization_.p(2,j)/dx2+ (discretization_.p(1,j-1)  + discretization_.p(1,j+1))/dy2-discretization_.rhs(1,j))/div_x;
            discretization_.setP(1,j, (1-omega) * discretization_.p(1,j)+ omega*  Gaussseidelterm1j);
                    

            
             for (int i = 2; i < i_max-2 ; i++)
            {
            //Gaussseideltermj1= v*((discretization_.p(i-1,1)+discretization_.p(i+1,1))/dx2+ discretization_.p(i,2)/ dy2-discretization_.rhs(i,1))/div_y;      
            //discretization_.setP(i,1, (1-omega) * discretization_.p(i,1)+ omega*  Gaussseideltermj1);
            
            x_term=(discretization_.p(i-1,j) + discretization_.p(i+1,j))/deltax_quad;
            y_term= (discretization_.p(i,j-1) + discretization_.p(i,j+1)) / deltay_quad;
            // value of pij gets overwritten with the new approximation
            discretization_.setP(i,j,(1-omega) * discretization_.p(i,j) + omega * vorfaktor*( x_term + y_term  - discretization_.rhs(i,j))) ;
            
              //Gaussseideltermijend= v*((discretization_.p(i-1,jend)+discretization_.p(i+1,jend))/dx2+ discretization_.p(i,jend-1)/dy2-discretization_.rhs(i,jend))/div_y;
              //discretization_.setP(i,jend, (1-omega) * discretization_.p(i,jend)+ omega* Gaussseideltermijend);
           
            } 
            
              Gaussseideltermiendj= v*(discretization_.p(iend-1,j)/dx2+ (discretization_.p(iend,j-1)+discretization_.p(iend,j+1))/dy2-discretization_.rhs(iend,j))/div_x;     
              discretization_.setP(iend,j, (1-omega) * discretization_.p(iend,j)+ omega* Gaussseideltermiendj);
            
       

        } 
        safe++;   
       discretization_.updatedPressureBC();
    }while(residuum() > discretization_.getepsilon() && safe<2000);
    std::cout<< "Residuum " << residuum() << " Safe "<< safe <<std::endl;
    //setPressureBoundaries();
     //discretization_.updatedPressureBC();
    }


