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
    {//FieldVariable p=discretization_.p();

    double deltax_quad = discretization_.dx() * discretization_.dx();
    double deltay_quad = discretization_.dy() * discretization_.dy();
    double vorfaktor= deltax_quad * deltay_quad/ (2 * (deltay_quad + deltax_quad ));
   // FieldVariable p=discretization_.p();
    int safe=0;

    //auxilliarys for code readability
    double x_term;
    double y_term;
    double omega=discretization_.getOmega();

    double epsilonquad=discretization_.getepsilon() *discretization_.getepsilon() ;   
    double resterm;
    int Nnumber= (discretization_.nCells()[0]*discretization_.nCells()[1]);

    do
    { //FieldVariable p=discretization_.p(); 
       
       
        for (int j = discretization_.pJBegin()+1; j < discretization_.pJEnd()-1; j++)

        {

        for (int i = discretization_.pIBegin()+1; i < discretization_.pIEnd()-1; i++)
            {   
            //discretization_.setP(i,1, v*((discretization_.p(i-1,1)+discretization_.p(i+1,1))/dx2+       discretization_.p(i,2)/dy2-discretization_.rhs(i,1))/div_y);
        

            x_term= (discretization_.p(i-1,j) + discretization_.p(i+1,j))/deltax_quad;
            y_term= (discretization_.p(i,j-1) + discretization_.p(i,j+1)) / deltay_quad;  
              // value of pij gets overwritten with the new approximation
          
            discretization_.setP(i,j, vorfaktor*( x_term + y_term   - discretization_.rhs(i,j))) ;
         
         
               //discretization_.setP(i,jend, v*((discretization_.p(i-1,jend)+discretization_.p(i+1,jend))/dx2+ discretization_.p(i,jend-1)/dy2-discretization_.rhs(i,jend))/div_y);
         
            
            }  
                 
           
        }        
         safe++;
      
         
        // discretization_.setP(p);
      // Ohne den konvergiert er doppelt so schnell  
      discretization_.updatedPressureBC();
         resterm=(residuum()*residuum())/Nnumber;
         
         
    }while(resterm > epsilonquad  && safe<20000);
    std::cout<< "Residuum " << residuum() << " Safe "<< safe <<std::endl;
    //discretization_.updatedPressureBC();
        
    
    }
    
/*/    
    
    
    void GaussSeidel::calculateP()
        {//FieldVariable p=discretization_.p();

    double deltax_quad = discretization_.dx() * discretization_.dx();
    double deltay_quad = discretization_.dy() * discretization_.dy();
    double vorfaktor= deltax_quad * deltay_quad/ (2 * (deltay_quad + deltax_quad ));
   // FieldVariable p=discretization_.p();
    int safe=0;

    //auxilliarys for code readability
    double x_term;
    double y_term;
    double omega=discretization_.getOmega();

    double epsilonquad=discretization_.getepsilon() *discretization_.getepsilon() ;   
    double resterm;
    int Nnumber= (discretization_.nCells()[0]*discretization_.nCells()[1]);

    double resterm_glob;
    double resterm_loc;
    
    do
    {  
       //Black part
       //Die erste Schleife iteriert ueber alle geraden Zeilen und geraden spalten
       //Die 2. Schleife ueber alle ungeraden zeilen und ungeraden Spalten
        for (int j = discretization_.pJBegin()+1; j < discretization_.pJEnd()-1; j+=2)
        {
          for (int i = discretization_.pIBegin()+1; i < discretization_.pIEnd()-1; i+=2)
              {   
              //discretization_.setP(i,1, v*((discretization_.p(i-1,1)+discretization_.p(i+1,1))/dx2+ discretization_.p(i,2)/dy2-discretization_.rhs(i,1))/div_y);
              x_term= (discretization_.p(i-1,j) + discretization_.p(i+1,j))/deltax_quad;
              y_term= (discretization_.p(i,j-1) + discretization_.p(i,j+1)) / deltay_quad;  
              discretization_.setP(i,j, vorfaktor*( x_term + y_term   - discretization_.rhs(i,j))) ;
        
             }    
        }   
        for (int j = discretization_.pJBegin()+2; j < discretization_.pJEnd()-1; j+=2)
        {
          for (int i = discretization_.pIBegin()+2; i < discretization_.pIEnd()-1; i+=2)
              {   
              //discretization_.setP(i,1, v*((discretization_.p(i-1,1)+discretization_.p(i+1,1))/dx2+ discretization_.p(i,2)/dy2-discretization_.rhs(i,1))/div_y);
              x_term= (discretization_.p(i-1,j) + discretization_.p(i+1,j))/deltax_quad;
              y_term= (discretization_.p(i,j-1) + discretization_.p(i,j+1)) / deltay_quad;  
              discretization_.setP(i,j, vorfaktor*( x_term + y_term   - discretization_.rhs(i,j))) ;
              }  
                  
        }
    //  discretization_.setPressureBCParalell(); // halbe Iteration vorbei,update pressureboundaries
         //   discretization_.updatedPressureBC();

      // Red Part
      //Erste Schleife geht ueber alle ungeraden Spalten und geraden Zeilen
      //Zweite Schleife ueber alle ungeraden Zeilen und geraden Spalten
       for (int j = discretization_.pJBegin()+2; j < discretization_.pJEnd()-1; j+=2)
        {
          for (int i = discretization_.pIBegin()+1; i < discretization_.pIEnd()-1; i+=2)
              {   
              //discretization_.setP(i,1, v*((discretization_.p(i-1,1)+discretization_.p(i+1,1))/dx2+ discretization_.p(i,2)/dy2-discretization_.rhs(i,1))/div_y);
              x_term= (discretization_.p(i-1,j) + discretization_.p(i+1,j))/deltax_quad;
              y_term= (discretization_.p(i,j-1) + discretization_.p(i,j+1)) / deltay_quad;  
              discretization_.setP(i,j, vorfaktor*( x_term + y_term   - discretization_.rhs(i,j))) ;
        
             }    
        }   
        for (int j = discretization_.pJBegin()+1; j < discretization_.pJEnd()-1; j+=2)
        {
          for (int i = discretization_.pIBegin()+2; i < discretization_.pIEnd()-1; i+=2)
              {   
              //discretization_.setP(i,1, v*((discretization_.p(i-1,1)+discretization_.p(i+1,1))/dx2+ discretization_.p(i,2)/dy2-discretization_.rhs(i,1))/div_y);
              x_term= (discretization_.p(i-1,j) + discretization_.p(i+1,j))/deltax_quad;
              y_term= (discretization_.p(i,j-1) + discretization_.p(i,j+1)) / deltay_quad;  
              discretization_.setP(i,j, vorfaktor*( x_term + y_term   - discretization_.rhs(i,j))) ;
              }  
                  
        }
//      discretization_.setPressureBCParalell(); // halbe Iteration vorbei,update pressureboundaries
        discretization_.updatedPressureBC();
    


                
         safe++;
      
         
       
      resterm=(residuum()*residuum())/Nnumber;

    //resterm_loc=(residuum()*residuum())/Nnumber;
     
    //MPI_Allreduce(&resterm_loc,&resterm_glob,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    //resterm=resterm_glob;


    }while(resterm > epsilonquad  && safe<discretization_.getMaxIteration());
    std::cout<< "Residuum " << residuum() << " Safe "<< safe <<std::endl;
   // printf("Rank %d Residuum %f Safe %d \n",discretization_.getOwnRankNo(),residuum(),safe);
    //discretization_.updatedPressureBC();    
    }

    */
