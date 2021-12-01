#include "discretization.h"
#include "staggeredgrid.h"
#include "settings.h"
#include "fieldvariable.h"


//constructor
Discretization::Discretization(Settings settings, Partitioning partitioning):StaggeredGrid(settings, partitioning) //maybe replace settings nCell with partitioning nCells.
,F({partitioning.nCells()[0]+3 - partitioning.ownPartitionContainsRightBoundary(), partitioning.nCells()[1]+3}, {0,0.5}, {settings.physicalSize[0] / (1.0*partitioning.nCells()[0]), settings.physicalSize[1] / (1.0*partitioning.nCells()[1])}) // is actually smaller than this, but makes handling easier
,G({partitioning.nCells()[0]+3, partitioning.nCells()[1]+3 - partitioning.ownPartitionContainsTopBoundary()}, {0.5,0}, {settings.physicalSize[0] / (1.0*partitioning.nCells()[0]), settings.physicalSize[1] / (1.0*partitioning.nCells()[1])}) // is actually smaller than this, but makes handling easier
,rhs_({partitioning.nCells()[0]+3, partitioning.nCells()[1]+3}, {0.5,0.5}, {settings.physicalSize[0] / (1.0*partitioning.nCells()[0]), settings.physicalSize[1] / (1.0*partitioning.nCells()[1])}) // is actually smaller than this, but makes handling easier
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
   
    if (partitioning_.ownPartitionContainsRightBoundary() )
    {  for (int j = pJBegin(); j < pJEnd(); j++)
            {
            pressure(i_max-1,j)=pressure(i_max-2,j);
            }

    }else   
    {      std::vector<double> Buffer_send(lj,0);
           std::vector<double> Buffer_recv(lj,0);
           int rank=partitioning_.coordiantesToRank(self_i+1,self_j);  
        
        for (int j = pJBegin(); j < pJEnd(); j++)
            {
            Buffer_send[j-pJBegin()]=pressure(i_max-2,j);          
            }
        if (partitioning_.first())   
             {
             MPI_Send(Buffer_send.data(),lj,MPI_DOUBLE,rank,3,MPI_COMM_WORLD);
             MPI_Recv(Buffer_recv.data(),lj,MPI_DOUBLE,rank,4,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
             }else
            {
             MPI_Recv(Buffer_recv.data(),lj,MPI_DOUBLE,rank,4,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
             MPI_Send(Buffer_send.data(),lj,MPI_DOUBLE,rank,3,MPI_COMM_WORLD);
             }
                
        for (int j = pJBegin(); j < j_max; j++)
            {
             pressure(i_max-1,j)=Buffer_recv[j-pJBegin()];
   
            }
    }    
    if (partitioning_.ownPartitionContainsLeftBoundary()) 
     {
           for (int j = pJBegin(); j < pJEnd(); j++)
            {
             pressure(pIBegin(),j)=pressure(pIBegin()+1,j);
            }
     }else
    {   std::vector<double> Buffer_send(lj,0);
        std::vector<double> Buffer_recv(lj,0);
        int rank=partitioning_.coordiantesToRank(self_i-1,self_j);  
        
       // std::cout<<"notleft " << rank <<std::endl;

        for (int j = pJBegin(); j < pJEnd(); j++)
            {
            Buffer_send[j-pJBegin()]=pressure(pJBegin()+1,j);           
            }
         if (partitioning_.first())   
         {
           MPI_Send(Buffer_send.data(),lj,MPI_DOUBLE,rank,4,MPI_COMM_WORLD);
           MPI_Recv(Buffer_recv.data(),lj,MPI_DOUBLE,rank,3,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
         } else {
           MPI_Recv(Buffer_recv.data(),lj,MPI_DOUBLE,rank,3,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
           MPI_Send(Buffer_send.data(),lj,MPI_DOUBLE,rank,4,MPI_COMM_WORLD);
         } 
          
          for (int j = pJBegin(); j < pJEnd(); j++)
            {
            pressure(pIBegin(),j)=Buffer_recv[j-pJBegin()];
            } 
    } 
    
    
////////////////////////////////////////////////////////////////////////////////
      if (partitioning_.ownPartitionContainsTopBoundary())
       {              
        for (int i = pIBegin(); i < pIEnd(); i++)
            {
               pressure(i,j_max-1)=pressure(i,j_max-2);
            }  
        }else 
        { 

            std::vector<double> Buffer_send(li,0);
            std::vector<double> Buffer_recv(li,0);
            int rank=partitioning_.coordiantesToRank(self_i,self_j+1);     
            
            for (int i = pIBegin(); i < pIEnd(); i++)
                {
                    Buffer_send[i-pIBegin()]=pressure(i,j_max-2);
                }
            if (partitioning_.first()) 
            {
            MPI_Send(Buffer_send.data(),li,MPI_DOUBLE,rank,5,MPI_COMM_WORLD);
            MPI_Recv(Buffer_recv.data(),li,MPI_DOUBLE,rank,6,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);          
            }else
            {
            MPI_Recv(Buffer_recv.data(),li,MPI_DOUBLE,rank,6,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);
            MPI_Send(Buffer_send.data(),li,MPI_DOUBLE,rank,5,MPI_COMM_WORLD);           
            }
           
            for (int i = pIBegin(); i < pIEnd(); i++)
                {
                  pressure(i,j_max-1)=Buffer_recv[i-pIBegin()];
                }
        } 
   if (partitioning_.ownPartitionContainsBottomBoundary())
            { 
            for (int i = pIBegin(); i < pIEnd(); i++)
               {
               pressure(i,pJBegin())=pressure(i,pJBegin()+1);                            
               }
           }else
           {
            std::vector<double> Buffer_send(li,0);
            std::vector<double> Buffer_recv(li,0);
            int rank=partitioning_.coordiantesToRank(self_i,self_j-1);    
                   
            for (int i = pIBegin(); i < pIEnd(); i++)
                {
                    Buffer_send[i-pIBegin()]=pressure(i,pJBegin()+1);
                }
            if (partitioning_.first()) 
            {
              MPI_Send(Buffer_send.data(),li,MPI_DOUBLE,rank,6,MPI_COMM_WORLD);
              MPI_Recv(Buffer_recv.data(),li,MPI_DOUBLE,rank,5,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);           
            }else
            {
              MPI_Recv(Buffer_recv.data(),li,MPI_DOUBLE,rank,5,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);
              MPI_Send(Buffer_send.data(),li,MPI_DOUBLE,rank,6,MPI_COMM_WORLD);           
            }
            for (int i = pIBegin(); i < pIEnd(); i++)
                {
                  pressure(i,pJBegin())=Buffer_recv[i-pIBegin()];
                }
        } 
}  



void Discretization::setBorderVelocityParalell(std::array<double,2> top,std::array<double,2> left,std::array<double,2> right,std::array<double,2> bottom)
{   int self_i=partitioning_.nodeOffset()[0];
    int self_j=partitioning_.nodeOffset()[1];
   
   //u left right
   int i_u_max  = velocity_X.size()[0], j_u_max = velocity_X.size()[1];
   int i_v_max  = velocity_Y.size()[0], j_v_max = velocity_Y.size()[1];
   int i_u_begin  = uIBegin(), j_u_begin = uJBegin();
   int i_v_begin  = vIBegin(), j_v_begin = vJBegin();
   int lu=(j_u_max-1)-(uJBegin()+1);
   int lv=(j_v_max-1)-(vJBegin()+1);
      if (partitioning_.ownPartitionContainsRightBoundary() )
   {       for (int j = uJBegin(); j < j_u_max; j++)
            {
                   // velocity_X(i_u_begin,j)=left[0];
                    velocity_X(i_u_max-1,j)=right[0];
            }
                for (int j = vJBegin(); j < j_v_max; j++)
            {
               // velocity_Y(i_v_begin,j)=2*left[1]-velocity_Y(i_v_begin+1,j); 
            velocity_Y(i_v_max-1,j)=2*right[1]-velocity_Y(i_v_max-2,j);
            }
        }else{
        std::vector<double> Buffer_send_u(lu,0);
        std::vector<double> Buffer_send_v(lv,0);
        std::vector<double> Buffer_recv_u(lu,0);
        std::vector<double> Buffer_recv_v(lv,0);
        int rank=partitioning_.coordiantesToRank(self_i+1,self_j);  
        
              for (int j = uJBegin()+1; j < j_u_max-1; j++)
           {
            Buffer_send_u[j-(j_u_begin+1)]=velocity_X(i_u_max-3,j);
            // Buffer_send_u[lu+j-(uJBegin()+1)]=velocity_X(i_u_max-3,j);
           }
        for (int j = vJBegin()+1; j < j_v_max-1; j++)
           {
            Buffer_send_v[j-(vJBegin()+1)]= velocity_Y(i_v_max-2,j);
           }

        if (partitioning_.first()) 
        {
            MPI_Send(Buffer_send_u.data(),lu,MPI_DOUBLE,rank,7,MPI_COMM_WORLD);
            MPI_Recv(Buffer_recv_u.data(),lu,MPI_DOUBLE,rank,8,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
            MPI_Send(Buffer_send_v.data(),lv,MPI_DOUBLE,rank,9,MPI_COMM_WORLD);           
            MPI_Recv(Buffer_recv_v.data(),lv,MPI_DOUBLE,rank,10,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
        }else
        {
            MPI_Recv(Buffer_recv_u.data(),lu,MPI_DOUBLE,rank,8,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
            MPI_Send(Buffer_send_u.data(),lu,MPI_DOUBLE,rank,7,MPI_COMM_WORLD);
            MPI_Recv(Buffer_recv_v.data(),lv,MPI_DOUBLE,rank,10,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);     
            MPI_Send(Buffer_send_v.data(),lv,MPI_DOUBLE,rank,9,MPI_COMM_WORLD);            
        }
        
        for (int j = uJBegin()+1; j < j_u_max-1; j++)
        {
             velocity_X(i_u_max-1,j)=Buffer_recv_u[j-(uJBegin()+1)];
        } 
        for (int j = vJBegin()+1; j < j_v_max-1; j++)
         {
            velocity_Y(i_v_max-1,j)=Buffer_recv_v[j-(vJBegin()+1)];
        }
    } 
   // }else if (!partitioning_.ownPartitionContainsLeftBoundary() && partitioning_.ownPartitionContainsRightBoundary())
        
  if (partitioning_.ownPartitionContainsLeftBoundary())
     {       for (int j = uJBegin(); j < j_u_max; j++)
            {
                    velocity_X(i_u_begin,j)=left[0];
                    //velocity_X(i_u_max-1,j)=right[0];
            }
                for (int j = vJBegin(); j < j_v_max; j++)
            {
                velocity_Y(i_v_begin,j)=2*left[1]-velocity_Y(i_v_begin+1,j); 
                //velocity_Y(i_v_max-1,j)=2*right[1]-velocity_Y(i_v_max-2,j);
            }
     }else{
         std::vector<double> Buffer_send_u(lu,0);
        std::vector<double> Buffer_send_v(lv,0);
        std::vector<double> Buffer_recv_u(lu,0);
        std::vector<double> Buffer_recv_v(lv,0);
        int rank=partitioning_.coordiantesToRank(self_i-1,self_j);  
        
        
        for (int j = uJBegin()+1; j < j_u_max-1; j++)
           {
            Buffer_send_u[j-(uJBegin()+1)]=velocity_X(i_u_begin+2,j);
            }
        for (int j = vJBegin()+1; j < j_v_max-1; j++)
           {
             Buffer_send_v[j-(vJBegin()+1)]=velocity_Y(i_v_begin+1,j);
           }
        if (partitioning_.first()) 
        {
           MPI_Send(Buffer_send_u.data(),lu,MPI_DOUBLE,rank,8,MPI_COMM_WORLD);
           MPI_Recv(Buffer_recv_u.data(),lu,MPI_DOUBLE,rank,7,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
           MPI_Send(Buffer_send_v.data(),lv,MPI_DOUBLE,rank,10,MPI_COMM_WORLD);           
           MPI_Recv(Buffer_recv_v.data(),lv,MPI_DOUBLE,rank,9,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
        }else
        {
           MPI_Recv(Buffer_recv_u.data(),lu,MPI_DOUBLE,rank,7,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
           MPI_Send(Buffer_send_u.data(),lu,MPI_DOUBLE,rank,8,MPI_COMM_WORLD);
           MPI_Recv(Buffer_recv_v.data(),lv,MPI_DOUBLE,rank,9,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
           MPI_Send(Buffer_send_v.data(),lv,MPI_DOUBLE,rank,10,MPI_COMM_WORLD);            
        }
          
        for (int j = uJBegin()+1; j < j_u_max-1; j++)
        {
             velocity_X(i_u_begin,j)=Buffer_recv_u[j-(uJBegin()+1)];
            //  velocity_X(i_u_begin,j)=Buffer_recv_u[(lu-2)+j-(uJBegin()+1)];
        } 
        for (int j = vJBegin()+1; j < j_v_max-1; j++)
         {
            velocity_Y(i_v_begin,j)=Buffer_recv_v[j-(vJBegin()+1)];
        }
   
        
    } 
    // top-bottom -------------------------------------------------------------------------------------------
    lu=(uIEnd()-1)-(uIBegin()+1);
    lv=(vIEnd()-1)-(vIBegin()+1);

    if (partitioning_.ownPartitionContainsBottomBoundary())
    {
        for (int i = uIBegin() + 1; i < i_u_max-1; i++) // i starts at 1 and goes to i_u_max-1 so that the wall is the BC in corners
        {
            velocity_X(i,uJBegin())=2*bottom[0]-velocity_X(i,vJBegin()+1);
        }
        for (int i = vIBegin()+1; i < vIEnd()-1; i++) // i starts at 1 and goes to i_u_max-1 so that the wall is the BC in corners
        {
            velocity_Y(i,vJBegin())=bottom[1]; 
        }
    } else // send/recive at bottom
    {
        std::vector<double> Buffer_send_u_B(lu,0);
        std::vector<double> Buffer_send_v_B(lv,0);
        std::vector<double> Buffer_recv_u_B(lu,0);
        std::vector<double> Buffer_recv_v_B(lv,0);
        int rank_B=partitioning_.coordiantesToRank(self_i,self_j-1);

        for (int i = vIBegin()+1; i < vIEnd()-1; i++)
        {
            Buffer_send_u_B[i-(vIBegin()+1)] = velocity_X(i, vIBegin()+1);
            Buffer_send_v_B[i-(vIBegin()+1)] = velocity_Y(i, vJBegin()+2);
        }
        
        if (partitioning_.first()) 
        {
          MPI_Send(Buffer_send_u_B.data(),lu,MPI_DOUBLE,rank_B,1,MPI_COMM_WORLD);
          MPI_Send(Buffer_send_v_B.data(),lv,MPI_DOUBLE,rank_B,2,MPI_COMM_WORLD);           
          MPI_Recv(Buffer_recv_u_B.data(),lu,MPI_DOUBLE,rank_B,1,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
          MPI_Recv(Buffer_recv_v_B.data(),lv,MPI_DOUBLE,rank_B,2,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
        }else
        {
          MPI_Recv(Buffer_recv_u_B.data(),lu,MPI_DOUBLE,rank_B,1,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
          MPI_Recv(Buffer_recv_v_B.data(),lv,MPI_DOUBLE,rank_B,2,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
          MPI_Send(Buffer_send_u_B.data(),lu,MPI_DOUBLE,rank_B,1,MPI_COMM_WORLD);
          MPI_Send(Buffer_send_v_B.data(),lv,MPI_DOUBLE,rank_B,2,MPI_COMM_WORLD);            
        }

        for (int i = vIBegin()+1; i < vIEnd()-1; i++)
        {
            velocity_X(i, vJBegin()) = Buffer_send_u_B[i-(vIBegin()+1)];
            velocity_Y(i, vJBegin()) = Buffer_send_v_B[i-(vIBegin()+1)];
        }
    }
   
    if (partitioning_.ownPartitionContainsTopBoundary())
    {
        for (int i = uIBegin() + 1; i < i_u_max-1; i++) // i starts at 1 and goes to i_u_max-1 so that the wall is the BC in corners
        {
            velocity_X(i,j_u_max-1)=2*top[0]-velocity_X(i,j_u_max-2);
        }
        for (int i = vIBegin()+1; i < vIEnd()-1; i++) // i starts at 1 and goes to i_u_max-1 so that the wall is the BC in corners
        {
            velocity_Y(i,vJEnd()-1)=top[1]; 
        }
    } else // send/recive top
    {
        std::vector<double> Buffer_send_u_T(lu,0);
        std::vector<double> Buffer_send_v_T(lv,0);
        std::vector<double> Buffer_recv_u_T(lu,0);
        std::vector<double> Buffer_recv_v_T(lv,0);
        int rank_T=partitioning_.coordiantesToRank(self_i,self_j+1);

        for (int i = vIBegin()+1; i < vIEnd()-1; i++)
        {
            Buffer_send_u_T[i-(vIBegin()+1)] = velocity_X(i, vJEnd()-2);
            Buffer_send_v_T[i-(vIBegin()+1)] = velocity_Y(i, vJEnd()-3);
        }

        if (partitioning_.first()) 
        {
          MPI_Send(Buffer_send_u_T.data(),lu,MPI_DOUBLE,rank_T,1,MPI_COMM_WORLD);
          MPI_Send(Buffer_send_v_T.data(),lv,MPI_DOUBLE,rank_T,2,MPI_COMM_WORLD);           
          MPI_Recv(Buffer_recv_u_T.data(),lu,MPI_DOUBLE,rank_T,1,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
          MPI_Recv(Buffer_recv_v_T.data(),lv,MPI_DOUBLE,rank_T,2,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
        }else
        {
          MPI_Recv(Buffer_recv_u_T.data(),lu,MPI_DOUBLE,rank_T,1,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
          MPI_Recv(Buffer_recv_v_T.data(),lv,MPI_DOUBLE,rank_T,2,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
          MPI_Send(Buffer_send_u_T.data(),lu,MPI_DOUBLE,rank_T,1,MPI_COMM_WORLD);
          MPI_Send(Buffer_send_v_T.data(),lv,MPI_DOUBLE,rank_T,2,MPI_COMM_WORLD);            
        }

        for (int i = vIBegin()+1; i < vIEnd()-1; i++)
        {
            velocity_X(i, uJEnd()-1) = Buffer_recv_u_T[i-(vIBegin()+1)]; //hier angepasst
            velocity_Y(i, vJEnd()-1) = Buffer_recv_v_T[i-(vIBegin()+1)];
         //std::cout<< (vIEnd()-1) -(vIBegin()+1) << " "  << (uIEnd()-1)- (uIBegin()+1) <<std::endl;        
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

