#include "discretization.h"
#include "staggeredgrid.h"
#include "settings.h"
#include "fieldvariable.h"


//constructor
Discretization::Discretization(Settings settings, Partitioning partitioning):StaggeredGrid(settings, partitioning) //maybe replace settings nCell with partitioning nCells.
,F({partitioning.nCells()[0]+3 - partitioning.ownPartitionContainsRightBoundary(), partitioning.nCells()[1]+3}, {1,1.5}, {settings.physicalSize[0] / (1.0*partitioning.nCellsGlobal()[0]), settings.physicalSize[1] / (1.0*partitioning.nCellsGlobal()[1])}) // is actually smaller than this, but makes handling easier
,G({partitioning.nCells()[0]+3, partitioning.nCells()[1]+3 - partitioning.ownPartitionContainsTopBoundary()}, {1.5,1}, {settings.physicalSize[0] / (1.0*partitioning.nCellsGlobal()[0]), settings.physicalSize[1] / (1.0*partitioning.nCellsGlobal()[1])}) // is actually smaller than this, but makes handling easier
,rhs_({partitioning.nCells()[0]+3, partitioning.nCells()[1]+3}, {1.5,1.5}, {settings.physicalSize[0] / (1.0*partitioning.nCellsGlobal()[0]), settings.physicalSize[1] / (1.0*partitioning.nCellsGlobal()[1])}) // is actually smaller than this, but makes handling easier
// ,partitioning_(partitioning)
{

}


// setBorderVelocity(settings.dirichletBcTop,settings.dirichletBcLeft,settings.dirichletBcRight,settings.dirichletBcBottom)
void Discretization::updateDeltaT()
{
    double time_limit_diffusion = (delta_x*delta_x*delta_y*delta_y)/((delta_x*delta_x)+(delta_y*delta_y))*settings_.re/2;
    double time_limit           = min3(time_limit_diffusion, delta_x/velocity_X.absmax(), delta_y/velocity_Y.absmax()) * settings_.tau;
    double deltat_loc=min2(time_limit, settings_.maximumDt);
  //   std::cout<< "Rank this  "<< partitioning_.ownRankNo()<< " " << deltat_loc <<std::endl;
    double deltat_glob=0; 
    MPI_Allreduce(&deltat_loc,&deltat_glob,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD); //Allreduce as every rank needs the global deltat
    deltat=deltat_glob;

   // std::cout << deltat <<std::endl;
}
void Discretization::setDeltaT(double delta_t)
{
   deltat=delta_t;
}

//destructor
Discretization::~Discretization()
{
}


// second order numerical derivation terms
double Discretization::computeDuDx2(int i, int j) const
{
 return (u(i+1,j) - 2 * u(i,j) + u(i-1,j))/(delta_x*delta_x);
}
double Discretization::computeDvDy2(int i, int j) const
{
 return (v(i,j+1) - 2 * v(i,j) + v(i,j-1))/ (delta_y*delta_y);
}

double Discretization::computeDuDy2(int i, int j) const
{
 return (u(i,j+1) - 2 * u(i,j) + u(i,j-1))/(delta_y*delta_y);
}
double Discretization::computeDvDx2(int i, int j) const
{
 return (v(i+1,j) - 2 * v(i,j) + v(i-1,j))/ (delta_x* delta_x);       
}


// firstt order numerical derivation terms
double Discretization::computeDuDx(int i, int j) const
{ 
  return (u(i,j)-u(i-1,j))/delta_x;
} 
double Discretization::computeDvDy(int i, int j) const
{ 
  return (v(i,j)-v(i,j-1))/delta_y;
} 

double Discretization::computeDpDx(int i, int j) const
{
  return (p(i+1, j)-p(i,j))/delta_x;
}

double Discretization::computeDpDy(int i, int j) const
{
  return (p(i, j+1)-p(i,j))/delta_y;
}


void Discretization::updateVelocity()
{   // Does not go over Boundary
    for (int j = uJBegin()+1; j < uJEnd()-1; j++)
    {
        for (int i = uIBegin()+1; i < uIEnd()-1; i++)
        {
            velocity_X(i, j)=F(i, j)-deltat*computeDpDx(i, j);
        }
        
    }
    for (int j = vJBegin()+1; j < vJEnd()-1; j++)
    {
        for (int i = vIBegin()+1; i < vIEnd()-1; i++)
        {
            velocity_Y(i, j)=G(i, j)-deltat*computeDpDy(i, j);
        }
        
    }
    
}

void Discretization::updateBoundaryFG()
{
    for (int j = uJBegin(); j < uJEnd(); j++)
   {
        F(uIBegin(),j)=u(uIBegin(),j);
        F(uIEnd()-1,j)=u(uIEnd()-1,j);
   }

   for (int i = vIBegin(); i < vIEnd(); i++)
   {
        G(i,vJBegin())=v(i,vJBegin());
        G(i,vJEnd()-1)=v(i,vJEnd()-1);

   }
}
        

double Discretization::min2(double value1, double value2) const
{
    if (value1<=value2)
    {
        return value1;
    }else 
    {
        return value2;
    }
}

double Discretization::min3(double value1, double value2, double value3) const
{
    if (value1<=value2 && value1<=value3)
    {
        return value1;
    }else if (value2<=value1 && value2<=value3)
    {
        return value2;
    } else
    {
        return value3;
    }
}

double Discretization::f(int i, int j) const
{
    return F(i,j);
}
double Discretization::g(int i, int j) const
{
    return G(i,j);
}
double Discretization::rhs(int i, int j) const//const statment removed because of error
{
    return rhs_(i, j);
}
double Discretization::getOmega() const
{
   return settings_.omega; 
} 

double Discretization::getDeltaT() const
{
   return deltat;
}


void Discretization::setRHS(int i, int j, double value)
{
    rhs_(i, j)=value;
}


void Discretization::setPressureBCParalell()
{   int self_i=partitioning_.nodeOffset()[0];
    int self_j=partitioning_.nodeOffset()[1];
        
    int i_max = pIEnd(), j_max = pJEnd();
    
    int lj=pJEnd()-pJBegin();
    int li=pIEnd()-pIBegin();
   
    if (partitioning_.ownPartitionContainsRightBoundary() &&  partitioning_.ownPartitionContainsLeftBoundary() )
    {  for (int j = pJBegin(); j < pJEnd(); j++)
            {
            pressure(i_max-1,j)=pressure(i_max-2,j);
            pressure(pIBegin(),j)=pressure(pIBegin()+1,j);
            }
    }
    else if(!partitioning_.ownPartitionContainsRightBoundary() &&  partitioning_.ownPartitionContainsLeftBoundary() )    
    {      std::vector<double> Buffer_send_R(lj,0);
           std::vector<double> Buffer_recv_R(lj,0);
           int rank=partitioning_.coordiantesToRank(self_i+1,self_j);  
        
        for (int j = pJBegin(); j < pJEnd(); j++)
            {
            Buffer_send_R[j-pJBegin()]=pressure(i_max-2,j);
            
             pressure(pIBegin(),j)=pressure(pIBegin()+1,j);
        
            //Buffer_send[j-pJBegin()]=100;
            //pressure(i_max-2,j)=400;          
            }
        if (partitioning_.first())   
             {
             MPI_Send(Buffer_send_R.data(),lj,MPI_DOUBLE,rank,3,MPI_COMM_WORLD);
             MPI_Recv(Buffer_recv_R.data(),lj,MPI_DOUBLE,rank,4,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
             }else
            {
             MPI_Recv(Buffer_recv_R.data(),lj,MPI_DOUBLE,rank,4,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
             MPI_Send(Buffer_send_R.data(),lj,MPI_DOUBLE,rank,3,MPI_COMM_WORLD);
             }
                
        for (int j = pJBegin(); j < j_max; j++)
            {
             pressure(i_max-1,j)=Buffer_recv_R[j-pJBegin()];
             //pressure(i_max-1,j)=300;
            }
    }    
    else if(!partitioning_.ownPartitionContainsLeftBoundary() && partitioning_.ownPartitionContainsRightBoundary()) 
     {
        std::vector<double> Buffer_send_L(lj,0);
        std::vector<double> Buffer_recv_L(lj,0);
        
        int rank=partitioning_.coordiantesToRank(self_i-1,self_j);  
        
       // std::cout<<"notleft " << rank <<std::endl;

        for (int j = pJBegin(); j < pJEnd(); j++)
            {
            Buffer_send_L[j-pJBegin()]=pressure(pJBegin()+1,j);
            pressure(i_max-1,j)=pressure(i_max-2,j);
            }

         if (partitioning_.first())   
         {
           MPI_Send(Buffer_send_L.data(),lj,MPI_DOUBLE,rank,4,MPI_COMM_WORLD);
           MPI_Recv(Buffer_recv_L.data(),lj,MPI_DOUBLE,rank,3,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
         } else {
           MPI_Recv(Buffer_recv_L.data(),lj,MPI_DOUBLE,rank,3,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
           MPI_Send(Buffer_send_L.data(),lj,MPI_DOUBLE,rank,4,MPI_COMM_WORLD);
         } 
          
          for (int j = pJBegin(); j < pJEnd(); j++)
            {
            pressure(pIBegin(),j)=Buffer_recv_L[j-pJBegin()];
            //pressure(pIBegin(),j)=300;
            } 

    } else if(!partitioning_.ownPartitionContainsLeftBoundary() && !partitioning_.ownPartitionContainsRightBoundary()) 
     {
        std::vector<double> Buffer_send_L(lj,0);
        std::vector<double> Buffer_recv_L(lj,0); 
        std::vector<double> Buffer_send_R(lj,0);
        std::vector<double> Buffer_recv_R(lj,0);
       
        int rank_L=partitioning_.coordiantesToRank(self_i-1,self_j);  
        int rank_R=partitioning_.coordiantesToRank(self_i+1,self_j);  
        
       // std::cout<<"notleft " << rank <<std::endl;

        for (int j = pJBegin(); j < pJEnd(); j++)
            {
            Buffer_send_L[j-pJBegin()]=pressure(pJBegin()+1,j);
            Buffer_send_R[j-pJBegin()]=pressure(i_max-2,j);       
            }
         if (partitioning_.first())   
         {
           MPI_Send(Buffer_send_R.data(),lj,MPI_DOUBLE,rank_R,3,MPI_COMM_WORLD);
           MPI_Send(Buffer_send_L.data(),lj,MPI_DOUBLE,rank_L,4,MPI_COMM_WORLD);
           MPI_Recv(Buffer_recv_R.data(),lj,MPI_DOUBLE,rank_R,4,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
           MPI_Recv(Buffer_recv_L.data(),lj,MPI_DOUBLE,rank_L,3,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
             
         }else {
           MPI_Recv(Buffer_recv_R.data(),lj,MPI_DOUBLE,rank_R,4,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
           MPI_Recv(Buffer_recv_L.data(),lj,MPI_DOUBLE,rank_L,3,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
           MPI_Send(Buffer_send_R.data(),lj,MPI_DOUBLE,rank_R,3,MPI_COMM_WORLD);
           MPI_Send(Buffer_send_L.data(),lj,MPI_DOUBLE,rank_L,4,MPI_COMM_WORLD);
         } 

           for (int j = pJBegin(); j < pJEnd(); j++)
            {
            pressure(pIBegin(),j)=Buffer_recv_L[j-pJBegin()];
            pressure(i_max-1,j)=Buffer_recv_R[j-pJBegin()];

            //pressure(pIBegin(),j)=300;
            } 

    } 
    
    
////////////////////////////////////////////////////////////////////////////////
     if (partitioning_.ownPartitionContainsTopBoundary() && partitioning_.ownPartitionContainsBottomBoundary())
       {              
        for (int i = pIBegin(); i < pIEnd(); i++)
            {
               pressure(i,j_max-1)=pressure(i,j_max-2);
               pressure(i,pJBegin())=pressure(i,pJBegin()+1);                  
            }  

        }else  if(!partitioning_.ownPartitionContainsTopBoundary() && partitioning_.ownPartitionContainsBottomBoundary())
        {   std::vector<double> Buffer_send_T(li,0);
            std::vector<double> Buffer_recv_T(li,0);
            int rank=partitioning_.coordiantesToRank(self_i,self_j+1);     
            
            for (int i = pIBegin(); i < pIEnd(); i++)
                {
                    Buffer_send_T[i-pIBegin()]=pressure(i,j_max-2);
                    //pressure(i,j_max-2)=400;
                   pressure(i,pJBegin())=pressure(i,pJBegin()+1);                  
                   
                }
            if (partitioning_.first()) 
            {
            MPI_Send(Buffer_send_T.data(),li,MPI_DOUBLE,rank,5,MPI_COMM_WORLD);
            MPI_Recv(Buffer_recv_T.data(),li,MPI_DOUBLE,rank,6,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);          
            }else
            {
            MPI_Recv(Buffer_recv_T.data(),li,MPI_DOUBLE,rank,6,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);
            MPI_Send(Buffer_send_T.data(),li,MPI_DOUBLE,rank,5,MPI_COMM_WORLD);           
            }
           
            for (int i = pIBegin(); i < pIEnd(); i++)
                {
                  pressure(i,j_max-1)=Buffer_recv_T[i-pIBegin()];
                 // pressure(i,j_max-1)=300;
                }
        } 
     else if (!partitioning_.ownPartitionContainsBottomBoundary() && partitioning_.ownPartitionContainsTopBoundary())
            { 
            std::vector<double> Buffer_send_B(li,0);
            std::vector<double> Buffer_recv_B(li,0);
            int rank=partitioning_.coordiantesToRank(self_i,self_j-1);    
                   
            for (int i = pIBegin(); i < pIEnd(); i++)
                {
                    Buffer_send_B[i-pIBegin()]=pressure(i,pJBegin()+1);
                    pressure(i,j_max-1)=pressure(i,j_max-2);
                }
            if (partitioning_.first()) 
            {
              MPI_Send(Buffer_send_B.data(),li,MPI_DOUBLE,rank,6,MPI_COMM_WORLD);
              MPI_Recv(Buffer_recv_B.data(),li,MPI_DOUBLE,rank,5,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);           
            }else
            {
              MPI_Recv(Buffer_recv_B.data(),li,MPI_DOUBLE,rank,5,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);
              MPI_Send(Buffer_send_B.data(),li,MPI_DOUBLE,rank,6,MPI_COMM_WORLD);           
            }
            for (int i = pIBegin(); i < pIEnd(); i++)
                {
                  pressure(i,pJBegin())=Buffer_recv_B[i-pIBegin()];
                  
                }

        } else if(!partitioning_.ownPartitionContainsBottomBoundary() && !partitioning_.ownPartitionContainsTopBoundary())
            { 
            std::vector<double> Buffer_send_B(li,0);
            std::vector<double> Buffer_recv_B(li,0);
            std::vector<double> Buffer_send_T(li,0);
            std::vector<double> Buffer_recv_T(li,0);
            
            int rank_B=partitioning_.coordiantesToRank(self_i,self_j-1);    
            int rank_T=partitioning_.coordiantesToRank(self_i,self_j+1);    
       
            for (int i = pIBegin(); i < pIEnd(); i++)
                {
                    Buffer_send_B[i-pIBegin()]=pressure(i,pJBegin()+1);
                    Buffer_send_T[i-pIBegin()]=pressure(i,j_max-2);
              
                    //pressure(i,pJBegin()+1)=400;
                 }
            if (partitioning_.first()) 
            { MPI_Send(Buffer_send_T.data(),li,MPI_DOUBLE,rank_T,5,MPI_COMM_WORLD);
              MPI_Send(Buffer_send_B.data(),li,MPI_DOUBLE,rank_B,6,MPI_COMM_WORLD);
              MPI_Recv(Buffer_recv_T.data(),li,MPI_DOUBLE,rank_T,6,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);
              MPI_Recv(Buffer_recv_B.data(),li,MPI_DOUBLE,rank_B,5,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);           
            }else
            { MPI_Recv(Buffer_recv_T.data(),li,MPI_DOUBLE,rank_T,6,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);          
              MPI_Recv(Buffer_recv_B.data(),li,MPI_DOUBLE,rank_B,5,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);
              MPI_Send(Buffer_send_T.data(),li,MPI_DOUBLE,rank_T,5,MPI_COMM_WORLD);           
              MPI_Send(Buffer_send_B.data(),li,MPI_DOUBLE,rank_B,6,MPI_COMM_WORLD);           
            }
    
            for (int i = pIBegin(); i < pIEnd(); i++)
                {
                  pressure(i,pJBegin())=Buffer_recv_B[i-pIBegin()];
                  pressure(i,j_max-1)=Buffer_recv_T[i-pIBegin()]; 
                  //pressure(i,pJBegin())=300;
                }

        }   
}  



void Discretization::setBorderVelocityParalell(std::array<double,2> top,std::array<double,2> left,std::array<double,2> right,std::array<double,2> bottom)
{   int self_i=partitioning_.nodeOffset()[0];
    int self_j=partitioning_.nodeOffset()[1];
   
   //u left right
   int i_u_max  = uIEnd(), j_u_max = uJEnd();
   int i_v_max  = vIEnd(), j_v_max = vJEnd();
   int i_u_begin  = uIBegin(), j_u_begin = uJBegin();
   int i_v_begin  = vIBegin(), j_v_begin = vJBegin();
   int lu=(j_u_max)-(uJBegin());
   int lv=(j_v_max)-(vJBegin());

      if (partitioning_.ownPartitionContainsRightBoundary() && partitioning_.ownPartitionContainsLeftBoundary())
   {       for (int j = uJBegin(); j < j_u_max; j++)
            {
                // velocity_X(i_u_begin,j)=left[0];
                velocity_X(i_u_max-1,j)=right[0];
                velocity_X(i_u_begin,j)=left[0];

            }
            for (int j = vJBegin(); j < j_v_max; j++)
            {
             // velocity_Y(i_v_begin,j)=2*left[1]-velocity_Y(i_v_begin+1,j); 
                velocity_Y(i_v_max-1,j)=2*right[1]-velocity_Y(i_v_max-2,j);
                velocity_Y(i_v_begin,j)=2*left[1]-velocity_Y(i_v_begin+1,j); 
              
            }

        }else if  (!partitioning_.ownPartitionContainsRightBoundary() && partitioning_.ownPartitionContainsLeftBoundary())
        {
        std::vector<double> Buffer_send_u_R(lu,0);
        std::vector<double> Buffer_send_v_R(lv,0);
        std::vector<double> Buffer_recv_u_R(lu,0);
        std::vector<double> Buffer_recv_v_R(lv,0);
        int rank=partitioning_.coordiantesToRank(self_i+1,self_j);  
        
              for (int j = uJBegin(); j < j_u_max; j++)
           {
            Buffer_send_u_R[j-(j_u_begin)]=velocity_X(i_u_max-3,j);
            velocity_X(i_u_begin,j)=left[0];
           }
        for (int j = vJBegin(); j < j_v_max; j++)
           {
            Buffer_send_v_R[j-(vJBegin())]= velocity_Y(i_v_max-2,j);
            velocity_Y(i_v_begin,j)=2*left[1]-velocity_Y(i_v_begin+1,j); 
           }
     
     
        if (partitioning_.first()) 
        {
            MPI_Send(Buffer_send_u_R.data(),lu,MPI_DOUBLE,rank,7,MPI_COMM_WORLD);
            MPI_Send(Buffer_send_v_R.data(),lv,MPI_DOUBLE,rank,9,MPI_COMM_WORLD);            
            MPI_Recv(Buffer_recv_u_R.data(),lu,MPI_DOUBLE,rank,8,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
            MPI_Recv(Buffer_recv_v_R.data(),lv,MPI_DOUBLE,rank,10,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
        }else
        {
            MPI_Recv(Buffer_recv_u_R.data(),lu,MPI_DOUBLE,rank,8,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
            MPI_Recv(Buffer_recv_v_R.data(),lv,MPI_DOUBLE,rank,10,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);       
            MPI_Send(Buffer_send_u_R.data(),lu,MPI_DOUBLE,rank,7,MPI_COMM_WORLD);
            MPI_Send(Buffer_send_v_R.data(),lv,MPI_DOUBLE,rank,9,MPI_COMM_WORLD);            
        }
    
        for (int j = uJBegin(); j < j_u_max; j++)
        {
             velocity_X(i_u_max-1,j)=Buffer_recv_u_R[j-(uJBegin())];
            //velocity_X(i_u_max-1,j)=100;
        } 
        for (int j = vJBegin(); j < j_v_max; j++)
         {
            velocity_Y(i_v_max-1,j)=Buffer_recv_v_R[j-(vJBegin())];
            //velocity_Y(i_v_max-1,j)=100;
        }
    } 
   // }else if (!partitioning_.ownPartitionContainsLeftBoundary() && partitioning_.ownPartitionContainsRightBoundary())
        
  if (!partitioning_.ownPartitionContainsLeftBoundary() && partitioning_.ownPartitionContainsRightBoundary())
     {       
        std::vector<double> Buffer_send_u_L(lu,0);
        std::vector<double> Buffer_send_v_L(lv,0);
        std::vector<double> Buffer_recv_u_L(lu,0);
        std::vector<double> Buffer_recv_v_L(lv,0);
        int rank=partitioning_.coordiantesToRank(self_i-1,self_j);  
        
        
        for (int j = uJBegin(); j < j_u_max; j++)
           {
            Buffer_send_u_L[j-(uJBegin())]=velocity_X(i_u_begin+2,j);
             velocity_X(i_u_max-1,j)=right[0];
            //velocity_X(i_u_begin+2,j)=200;
            }
        for (int j = vJBegin(); j < j_v_max; j++)
           {
             Buffer_send_v_L[j-(vJBegin())]=velocity_Y(i_v_begin+1,j);
             velocity_Y(i_v_max-1,j)=2*right[1]-velocity_Y(i_v_max-2,j);
            
            //velocity_Y(i_v_begin+2,j)=200;
           }
        if (partitioning_.first()) 
        {
           MPI_Send(Buffer_send_u_L.data(),lu,MPI_DOUBLE,rank,8,MPI_COMM_WORLD);
           MPI_Send(Buffer_send_v_L.data(),lv,MPI_DOUBLE,rank,10,MPI_COMM_WORLD);           
           MPI_Recv(Buffer_recv_u_L.data(),lu,MPI_DOUBLE,rank,7,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);           
           MPI_Recv(Buffer_recv_v_L.data(),lv,MPI_DOUBLE,rank,9,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
        }else
        {
           MPI_Recv(Buffer_recv_u_L.data(),lu,MPI_DOUBLE,rank,7,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
           MPI_Recv(Buffer_recv_v_L.data(),lv,MPI_DOUBLE,rank,9,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);              
           MPI_Send(Buffer_send_u_L.data(),lu,MPI_DOUBLE,rank,8,MPI_COMM_WORLD);
           MPI_Send(Buffer_send_v_L.data(),lv,MPI_DOUBLE,rank,10,MPI_COMM_WORLD);            
        }
          
        for (int j = uJBegin(); j < j_u_max; j++)
        {
             velocity_X(i_u_begin,j)=Buffer_recv_u_L[j-(uJBegin())];
        } 
        for (int j = vJBegin(); j < j_v_max; j++)
         {
            velocity_Y(i_v_begin,j)=Buffer_recv_v_L[j-(vJBegin())];
         }
   
        
    }else if(!partitioning_.ownPartitionContainsLeftBoundary() && !partitioning_.ownPartitionContainsRightBoundary())
     {       
        std::vector<double> Buffer_send_u_L(lu,0);
        std::vector<double> Buffer_send_v_L(lv,0);
        std::vector<double> Buffer_recv_u_L(lu,0);
        std::vector<double> Buffer_recv_v_L(lv,0);
        std::vector<double> Buffer_send_u_R(lu,0);
        std::vector<double> Buffer_send_v_R(lv,0);
        std::vector<double> Buffer_recv_u_R(lu,0);
        std::vector<double> Buffer_recv_v_R(lv,0);
        int rank_L=partitioning_.coordiantesToRank(self_i-1,self_j);  
        int rank_R=partitioning_.coordiantesToRank(self_i+1,self_j);  
        
        for (int j = uJBegin(); j < j_u_max; j++)
           {
            Buffer_send_u_L[j-(uJBegin())]=velocity_X(i_u_begin+2,j);
            Buffer_send_u_R[j-(j_u_begin)]=velocity_X(i_u_max-3,j);
            }
        for (int j = vJBegin(); j < j_v_max; j++)
           {
             Buffer_send_v_L[j-(vJBegin())]=velocity_Y(i_v_begin+1,j);
             Buffer_send_v_R[j-(vJBegin())]= velocity_Y(i_v_max-2,j);
           }
      
          
        if (partitioning_.first()) 
        {
           MPI_Send(Buffer_send_u_L.data(),lu,MPI_DOUBLE,rank_L,8,MPI_COMM_WORLD);
           MPI_Send(Buffer_send_u_R.data(),lu,MPI_DOUBLE,rank_R,7,MPI_COMM_WORLD);
           MPI_Send(Buffer_send_v_L.data(),lv,MPI_DOUBLE,rank_L,10,MPI_COMM_WORLD);            
           MPI_Send(Buffer_send_v_R.data(),lv,MPI_DOUBLE,rank_R,9,MPI_COMM_WORLD);            
           MPI_Recv(Buffer_recv_u_L.data(),lu,MPI_DOUBLE,rank_L,7,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);           
           MPI_Recv(Buffer_recv_u_R.data(),lu,MPI_DOUBLE,rank_R,8,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
           MPI_Recv(Buffer_recv_v_L.data(),lv,MPI_DOUBLE,rank_L,9,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);      
           MPI_Recv(Buffer_recv_v_R.data(),lv,MPI_DOUBLE,rank_R,10,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);     
          }else
        {
           MPI_Recv(Buffer_recv_u_L.data(),lu,MPI_DOUBLE,rank_L,7,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);           
           MPI_Recv(Buffer_recv_u_R.data(),lu,MPI_DOUBLE,rank_R,8,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
           MPI_Recv(Buffer_recv_v_L.data(),lv,MPI_DOUBLE,rank_L,9,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
           MPI_Recv(Buffer_recv_v_R.data(),lv,MPI_DOUBLE,rank_R,10,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);     
           MPI_Send(Buffer_send_u_L.data(),lu,MPI_DOUBLE,rank_L,8,MPI_COMM_WORLD);
           MPI_Send(Buffer_send_u_R.data(),lu,MPI_DOUBLE,rank_R,7,MPI_COMM_WORLD);
           MPI_Send(Buffer_send_v_L.data(),lv,MPI_DOUBLE,rank_L,10,MPI_COMM_WORLD);            
           MPI_Send(Buffer_send_v_R.data(),lv,MPI_DOUBLE,rank_R,9,MPI_COMM_WORLD);            
         }
          
        for (int j = uJBegin(); j < j_u_max; j++)
        {
             velocity_X(i_u_begin,j)=Buffer_recv_u_L[j-(uJBegin())];
             velocity_X(i_u_max-1,j)=Buffer_recv_u_R[j-(uJBegin())];

        } 
        for (int j = vJBegin(); j < j_v_max; j++)
         {
            velocity_Y(i_v_begin,j)=Buffer_recv_v_L[j-(vJBegin())];
            velocity_Y(i_v_max-1,j)=Buffer_recv_v_R[j-(vJBegin())];
        }
  
     }

    // top-bottom -------------------------------------------------------------------------------------------
    lu=(uIEnd()-1)-(uIBegin()+1);
    lv=(vIEnd()-1)-(vIBegin()+1);

    if (partitioning_.ownPartitionContainsBottomBoundary() && partitioning_.ownPartitionContainsTopBoundary())
    {
        for (int i = uIBegin()+partitioning_.ownPartitionContainsLeftBoundary(); i < i_u_max-partitioning_.ownPartitionContainsRightBoundary(); i++) // i starts at 1 and goes to i_u_max-1 so that the wall is the BC in corners
        {
            velocity_X(i,uJBegin())=2*bottom[0]-velocity_X(i,vJBegin()+1);
            velocity_X(i,j_u_max-1)=2*top[0]-velocity_X(i,j_u_max-2);
        
        }
        for (int i = vIBegin()+partitioning_.ownPartitionContainsLeftBoundary(); i < vIEnd()-partitioning_.ownPartitionContainsRightBoundary(); i++) // i starts at 1 and goes to i_u_max-1 so that the wall is the BC in corners
        {
            velocity_Y(i,vJBegin())=bottom[1];
            velocity_Y(i,vJEnd()-1)=top[1]; 
        }


    } else if (!partitioning_.ownPartitionContainsBottomBoundary() && partitioning_.ownPartitionContainsTopBoundary())
    {
        std::vector<double> Buffer_send_u_B(lu,0);
        std::vector<double> Buffer_send_v_B(lv,0);
        std::vector<double> Buffer_recv_u_B(lu,0);
        std::vector<double> Buffer_recv_v_B(lv,0);
        int rank_B=partitioning_.coordiantesToRank(self_i,self_j-1);

        for (int i = uIBegin()+1; i < uIEnd()-1; i++)
        {
            Buffer_send_u_B[i-(uIBegin()+1)] = velocity_X(i, uJBegin()+1);
            velocity_X(i,j_u_max-1)=2*top[0]-velocity_X(i,j_u_max-2);
  
        }
        
        for (int i = vIBegin()+1; i < vIEnd()-1; i++)
        {
            Buffer_send_v_B[i-(vIBegin()+1)] = velocity_Y(i, vJBegin()+2);
            velocity_Y(i,vJEnd()-1)=top[1]; 
        }
        
        if (partitioning_.first()) 
        {
          MPI_Send(Buffer_send_u_B.data(),lu,MPI_DOUBLE,rank_B,1,MPI_COMM_WORLD);
          MPI_Send(Buffer_send_v_B.data(),lv,MPI_DOUBLE,rank_B,2,MPI_COMM_WORLD);           
          MPI_Recv(Buffer_recv_u_B.data(),lu,MPI_DOUBLE,rank_B,11,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
          MPI_Recv(Buffer_recv_v_B.data(),lv,MPI_DOUBLE,rank_B,12,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
        }else
        {
          MPI_Recv(Buffer_recv_u_B.data(),lu,MPI_DOUBLE,rank_B,11,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
          MPI_Recv(Buffer_recv_v_B.data(),lv,MPI_DOUBLE,rank_B,12,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
          MPI_Send(Buffer_send_u_B.data(),lu,MPI_DOUBLE,rank_B,1,MPI_COMM_WORLD);
          MPI_Send(Buffer_send_v_B.data(),lv,MPI_DOUBLE,rank_B,2,MPI_COMM_WORLD);            
        }

        for (int i = uIBegin()+1; i < uIEnd()-1; i++)
        {
            velocity_X(i, uJBegin()) = Buffer_recv_u_B[i-(uIBegin()+1)];
        }

        for (int i = vIBegin()+1; i < vIEnd()-1; i++)
        {
            velocity_Y(i, vJBegin()) = Buffer_recv_v_B[i-(vIBegin()+1)];
        }
    }
    else if (!partitioning_.ownPartitionContainsTopBoundary() && partitioning_.ownPartitionContainsBottomBoundary() )
        {
      
        std::vector<double> Buffer_send_u_T(lu,0);
        std::vector<double> Buffer_send_v_T(lv,0);
        std::vector<double> Buffer_recv_u_T(lu,0);
        std::vector<double> Buffer_recv_v_T(lv,0);
        int rank_T=partitioning_.coordiantesToRank(self_i,self_j+1);

        for (int i = uIBegin()+1; i < uIEnd()-1; i++)
        {
            Buffer_send_u_T[i-(uIBegin()+1)] = velocity_X(i, uJEnd()-2);
        }

        for (int i = vIBegin()+1; i < vIEnd()-1; i++)
        {
            Buffer_send_v_T[i-(vIBegin()+1)] = velocity_Y(i, vJEnd()-3);
           
        }

        if (partitioning_.first()) 
        {
          MPI_Send(Buffer_send_u_T.data(),lu,MPI_DOUBLE,rank_T,11,MPI_COMM_WORLD);
          MPI_Send(Buffer_send_v_T.data(),lv,MPI_DOUBLE,rank_T,12,MPI_COMM_WORLD);           
          MPI_Recv(Buffer_recv_u_T.data(),lu,MPI_DOUBLE,rank_T,1,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
          MPI_Recv(Buffer_recv_v_T.data(),lv,MPI_DOUBLE,rank_T,2,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
        }else
        {
          MPI_Recv(Buffer_recv_u_T.data(),lu,MPI_DOUBLE,rank_T,1,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
          MPI_Recv(Buffer_recv_v_T.data(),lv,MPI_DOUBLE,rank_T,2,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
          MPI_Send(Buffer_send_u_T.data(),lu,MPI_DOUBLE,rank_T,11,MPI_COMM_WORLD);
          MPI_Send(Buffer_send_v_T.data(),lv,MPI_DOUBLE,rank_T,12,MPI_COMM_WORLD);            
        }
        for (int i = uIBegin()+1; i < uIEnd()-1; i++)
        {
            velocity_X(i, uJEnd()-1) = Buffer_recv_u_T[i-(uIBegin()+1)]; 
           // velocity_X(i, uJEnd()-1)=100;
        }
        for (int i = vIBegin()+1; i < vIEnd()-1; i++)
        {
            velocity_Y(i, vJEnd()-1) = Buffer_recv_v_T[i-(vIBegin()+1)];
            //velocity_Y(i, vJEnd()-1) = 100;
       
        } 
         

           
    }
    
  else if (!partitioning_.ownPartitionContainsTopBoundary() && !partitioning_.ownPartitionContainsBottomBoundary() )
        {
        std::vector<double> Buffer_send_u_B(lu,0);
        std::vector<double> Buffer_send_v_B(lv,0);
        std::vector<double> Buffer_recv_u_B(lu,0);
        std::vector<double> Buffer_recv_v_B(lv,0);
        int rank_B=partitioning_.coordiantesToRank(self_i,self_j-1);


        std::vector<double> Buffer_send_u_T(lu,0);
        std::vector<double> Buffer_send_v_T(lv,0);
        std::vector<double> Buffer_recv_u_T(lu,0);
        std::vector<double> Buffer_recv_v_T(lv,0);
        int rank_T=partitioning_.coordiantesToRank(self_i,self_j+1);

        for (int i = uIBegin()+1; i < uIEnd()-1; i++)
        {
            Buffer_send_u_T[i-(uIBegin()+1)] = velocity_X(i, uJEnd()-2);
            Buffer_send_u_B[i-(uIBegin()+1)] = velocity_X(i, uJBegin()+1);
            }

        for (int i = vIBegin()+1; i < vIEnd()-1; i++)
        {
            Buffer_send_v_T[i-(vIBegin()+1)] = velocity_Y(i, vJEnd()-3);
            Buffer_send_v_B[i-(vIBegin()+1)] = velocity_Y(i, vJBegin()+2);           
        }
   
       if (partitioning_.first()) 
        {
          MPI_Send(Buffer_send_u_B.data(),lu,MPI_DOUBLE,rank_B,1,MPI_COMM_WORLD);
          MPI_Send(Buffer_send_v_B.data(),lv,MPI_DOUBLE,rank_B,2,MPI_COMM_WORLD);           
   
          MPI_Send(Buffer_send_u_T.data(),lu,MPI_DOUBLE,rank_T,11,MPI_COMM_WORLD);
          MPI_Send(Buffer_send_v_T.data(),lv,MPI_DOUBLE,rank_T,12,MPI_COMM_WORLD);           
          
          MPI_Recv(Buffer_recv_u_B.data(),lu,MPI_DOUBLE,rank_B,11,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
          MPI_Recv(Buffer_recv_v_B.data(),lv,MPI_DOUBLE,rank_B,12,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    

          MPI_Recv(Buffer_recv_u_T.data(),lu,MPI_DOUBLE,rank_T,1,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
          MPI_Recv(Buffer_recv_v_T.data(),lv,MPI_DOUBLE,rank_T,2,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
        }else
        {
          MPI_Recv(Buffer_recv_u_B.data(),lu,MPI_DOUBLE,rank_B,11,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
          MPI_Recv(Buffer_recv_v_B.data(),lv,MPI_DOUBLE,rank_B,12,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    

          MPI_Recv(Buffer_recv_u_T.data(),lu,MPI_DOUBLE,rank_T,1,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
          MPI_Recv(Buffer_recv_v_T.data(),lv,MPI_DOUBLE,rank_T,2,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
       
          MPI_Send(Buffer_send_u_B.data(),lu,MPI_DOUBLE,rank_B,1,MPI_COMM_WORLD);
          MPI_Send(Buffer_send_v_B.data(),lv,MPI_DOUBLE,rank_B,2,MPI_COMM_WORLD);           
   
          MPI_Send(Buffer_send_u_T.data(),lu,MPI_DOUBLE,rank_T,11,MPI_COMM_WORLD);
          MPI_Send(Buffer_send_v_T.data(),lv,MPI_DOUBLE,rank_T,12,MPI_COMM_WORLD);           
        }   
   
   
        for (int i = uIBegin()+1; i < uIEnd()-1; i++)
        {
            velocity_X(i, uJEnd()-1) = Buffer_recv_u_T[i-(uIBegin()+1)]; 
            velocity_X(i, uJBegin()) = Buffer_recv_u_B[i-(uIBegin()+1)];
          
                 }
        for (int i = vIBegin()+1; i < vIEnd()-1; i++)
        {
            velocity_Y(i, vJEnd()-1) = Buffer_recv_v_T[i-(vIBegin()+1)];
            velocity_Y(i, vJBegin()) = Buffer_recv_v_B[i-(vIBegin()+1)];
       
        } 
        
          
        
        }
}      
   
   
   
   
   
   





void Discretization::updateBoundaryFGParalell()
{
    if (partitioning_.ownPartitionContainsLeftBoundary())
    {
        for (int j = uJBegin(); j < uJEnd(); j++)
        {
            F(uIBegin(),j)=u(uIBegin(),j);
        }
    }
    if (partitioning_.ownPartitionContainsRightBoundary())
    {
        for (int j = uJBegin(); j < uJEnd(); j++)
        {   
            F(uIEnd()-1,j)=u(uIEnd()-1,j);
        }   
    }

    if (partitioning_.ownPartitionContainsBottomBoundary())
    {
        for (int i = vIBegin(); i < vIEnd(); i++)
        {
            G(i,vJBegin())=v(i,vJBegin());
        }   
    }
    if (partitioning_.ownPartitionContainsTopBoundary())
    {
        for (int i = vIBegin(); i < vIEnd(); i++)
        {
            G(i,vJEnd()-1)=v(i,vJEnd()-1);
        }   
    }
   
}






int Discretization::getOwnRankNo()
{
    return partitioning_.ownRankNo();
}

