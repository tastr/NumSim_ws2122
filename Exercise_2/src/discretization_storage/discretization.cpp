#include "discretization.h"
#include "staggeredgrid.h"
#include "settings.h"
#include "fieldvariable.h"


//constructor
Discretization::Discretization(Settings settings, Partitioning partitioning):StaggeredGrid(settings) //maybe replace settings nCell with partitioning nCells.
,F({settings.nCells[0]+1, settings.nCells[1]+2}, {0,0.5}, {settings.physicalSize[0] / (1.0*settings.nCells[0]), settings.physicalSize[1] / (1.0*settings.nCells[1])}) // is actually smaller than this, but makes handling easier
,G({settings.nCells[0]+2, settings.nCells[1]+1}, {0.5,0}, {settings.physicalSize[0] / (1.0*settings.nCells[0]), settings.physicalSize[1] / (1.0*settings.nCells[1])}) // is actually smaller than this, but makes handling easier
,rhs_({settings.nCells[0]+2, settings.nCells[1]+2}, {0.5,0.5}, {settings.physicalSize[0] / (1.0*settings.nCells[0]), settings.physicalSize[1] / (1.0*settings.nCells[1])}) // is actually smaller than this, but makes handling easier
,partitioning_(partitioning)
{

}


// setBorderVelocity(settings.dirichletBcTop,settings.dirichletBcLeft,settings.dirichletBcRight,settings.dirichletBcBottom)
void Discretization::updateDeltaT()
{
    double time_limit_diffusion = (delta_x*delta_x*delta_y*delta_y)/((delta_x*delta_x)+(delta_y*delta_y))*settings_.re/2;
    double time_limit           = min3(time_limit_diffusion, delta_x/velocity_X.absmax(), delta_y/velocity_Y.absmax()) * settings_.tau;
    double deltat_loc=min2(time_limit, settings_.maximumDt);
    double deltat_glob; 
    MPI_Reduce(&deltat_loc,&deltat_glob,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
    deltat=deltat_glob;
}

//destructor
Discretization::~Discretization()
{
}

 double Discretization::computeDu2Dx(int i, int j) const {return 0;}   

 double Discretization::computeDv2Dy(int i, int j) const {return 0;}
 
 double Discretization::computeDuvDx(int i, int j) const {return 0;}
 
 double Discretization::computeDuvDy(int i, int j) const {return 0;}
 



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


// first order numerical derivation terms
double Discretization::computeDuDx(int i, int j) const
{ 
  return (u(i,j)-u(i-1,j))/delta_x;
} 
double Discretization::computeDvDy(int i, int j) const
{ 
  return (v(i,j)-v(i,j-1))/ delta_y;
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
        F(0,j)=u(0,j);
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


void Discretization::calculation()
{
     //this should never be called, as it is a virtual function
    assert(false);

   
}

void Discretization::setRHS(int i, int j, double value)
{
    rhs_(i, j)=value;
}


void Discretization::setPressureBCParalell()
{   int self_i=partitioning_.nodeOffset()[0];
    int self_j=partitioning_.nodeOffset()[1];
        
    int i_max = pressure.size()[0], j_max = pressure.size()[1];
    if (partitioning_.ownPartitionContainsLeftBoundary() && partitioning_.ownPartitionContainsRightBoundary())
        {
            for (int j = 0; j < j_max; j++)
            {
            pressure(0,j)=pressure(1,j);
            pressure(i_max-1,j)=pressure(i_max-2,j);
            }
    }else if (partitioning_.ownPartitionContainsLeftBoundary() && !partitioning_.ownPartitionContainsRightBoundary())
        {  std::vector<double> Buffer_send(j_max,0);
        std::vector<double> Buffer_recv(j_max,0);
        int rank=partitioning_.coordiantesToRank(self_i+1,self_j);  
        for (int j = 0; j < j_max; j++)
            {
            pressure(0,j)=pressure(1,j);
            Buffer_send[j]=pressure(i_max-2,j);
            }
        MPI_Send(Buffer_send.data(),j_max,MPI_DOUBLE,rank,1,MPI_COMM_WORLD);
        MPI_Recv(Buffer_recv.data(),j_max,MPI_DOUBLE,rank,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
        for (int j = 0; j < j_max; j++)
            {
            pressure(i_max-1,j)=Buffer_recv[j];
            } 
        }
    else if (!partitioning_.ownPartitionContainsLeftBoundary() && partitioning_.ownPartitionContainsRightBoundary())
    {   std::vector<double> Buffer_send(j_max,0);
        std::vector<double> Buffer_recv(j_max,0);
        int rank=partitioning_.coordiantesToRank(self_i+1,self_j);  
        for (int j = 0; j < j_max; j++)
            {
            pressure(i_max-1,j)=pressure(i_max-2,j);
            Buffer_send[j]=pressure(1,j);
            }
         if (partitioning_.ownRankNo()%2==0)   
         {
           MPI_Send(Buffer_send.data(),j_max,MPI_DOUBLE,rank,1,MPI_COMM_WORLD);
           MPI_Recv(Buffer_recv.data(),j_max,MPI_DOUBLE,rank,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
         } else
         {
          MPI_Recv(Buffer_recv.data(),j_max,MPI_DOUBLE,rank,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
          MPI_Send(Buffer_send.data(),j_max,MPI_DOUBLE,rank,1,MPI_COMM_WORLD);
         } 
          
          for (int j = 0; j < j_max; j++)
            {
            pressure(0,j)=Buffer_recv[j];
            } 
    }else
        {
        std::vector<double> Buffer_send_l(j_max,0);
        std::vector<double> Buffer_send_r(j_max,0);
        std::vector<double> Buffer_recv_l(j_max,0);
        std::vector<double> Buffer_recv_r(j_max,0);
        int rank_r=partitioning_.coordiantesToRank(self_i+1,self_j);  
        int rank_l=partitioning_.coordiantesToRank(self_i-1,self_j);  
        for (int j = 0; j < j_max; j++)
        {
            Buffer_send_l[j]=pressure(1,j);
            Buffer_send_r[j]=pressure(i_max-2,j);
        }
        if (partitioning_.ownRankNo()%2==0) 
        {
          MPI_Send(Buffer_send_r.data(),j_max,MPI_DOUBLE,rank_l,1,MPI_COMM_WORLD);
          MPI_Send(Buffer_send_r.data(),j_max,MPI_DOUBLE,rank_r,1,MPI_COMM_WORLD);           
          MPI_Recv(Buffer_recv_r.data(),j_max,MPI_DOUBLE,rank_r,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
          MPI_Recv(Buffer_recv_l.data(),j_max,MPI_DOUBLE,rank_l,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
       
        } else
        {
          MPI_Recv(Buffer_recv_r.data(),j_max,MPI_DOUBLE,rank_r,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
          MPI_Recv(Buffer_recv_l.data(),j_max,MPI_DOUBLE,rank_l,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
          MPI_Send(Buffer_send_r.data(),j_max,MPI_DOUBLE,rank_l,1,MPI_COMM_WORLD);
          MPI_Send(Buffer_send_r.data(),j_max,MPI_DOUBLE,rank_r,1,MPI_COMM_WORLD);               
        }
        for (int j = 0; j < j_max; j++)
         {
            pressure(0,j)=Buffer_recv_l[j];
            pressure(i_max-1,j)=Buffer_recv_r[j];
         }        
        }
////////////////////////////////////////////////////////////////////////////////
   if (partitioning_.ownPartitionContainsTopBoundary() && partitioning_.ownPartitionContainsBottomBoundary())
       for (int i = 0; i < i_max; i++)
        {
            pressure(i,0)=pressure(i,1);
            pressure(i,j_max-1)=pressure(i,j_max-2);
        }
    if (!partitioning_.ownPartitionContainsTopBoundary() && partitioning_.ownPartitionContainsBottomBoundary())
       { 
        std::vector<double> Buffer_send(i_max,0);
        std::vector<double> Buffer_recv(i_max,0);
        int rank=partitioning_.coordiantesToRank(self_i,self_j+1);    
        for (int i = 0; i < i_max; i++)
            {
                pressure(i,0)=pressure(i,1);
                Buffer_send[i]=pressure(i,j_max-2);
            }
        if (partitioning_.ownRankNo()%2==0) 
        {
         MPI_Send(Buffer_send.data(),i_max,MPI_DOUBLE,rank,1,MPI_COMM_WORLD);
         MPI_Recv(Buffer_recv.data(),i_max,MPI_DOUBLE,rank,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);          
        }else
        {
         MPI_Recv(Buffer_recv.data(),i_max,MPI_DOUBLE,rank,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);
         MPI_Send(Buffer_send.data(),i_max,MPI_DOUBLE,rank,1,MPI_COMM_WORLD);           
        }
        for (int i = 0; i < i_max; i++)
            {
                pressure(i,j_max-1)=Buffer_recv[i];
            }
        }    
        if (partitioning_.ownPartitionContainsTopBoundary() && !partitioning_.ownPartitionContainsBottomBoundary())
        {
        std::vector<double> Buffer_send(i_max,0);
        std::vector<double> Buffer_recv(i_max,0);
        int rank=partitioning_.coordiantesToRank(self_i,self_j-1);    
        for (int i = 0; i < i_max; i++)
            {
                pressure(i,j_max-1)=pressure(i,j_max-2);
                Buffer_send[i]=pressure(i,1);
            }
        if (partitioning_.ownRankNo()%2==0) 
        {
         MPI_Send(Buffer_send.data(),i_max,MPI_DOUBLE,rank,1,MPI_COMM_WORLD);
         MPI_Recv(Buffer_recv.data(),i_max,MPI_DOUBLE,rank,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);           
        }else
        {
         MPI_Recv(Buffer_recv.data(),i_max,MPI_DOUBLE,rank,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);
         MPI_Send(Buffer_send.data(),i_max,MPI_DOUBLE,rank,1,MPI_COMM_WORLD);           
        }
        for (int i = 0; i < i_max; i++)
            {
                pressure(i,j_max-1)=Buffer_recv[i];
            }
        }else
        {
        std::vector<double> Buffer_send_B(i_max,0);
        std::vector<double> Buffer_send_T(i_max,0);
        std::vector<double> Buffer_recv_B(i_max,0);
        std::vector<double> Buffer_recv_T(i_max,0);
        int rank_T=partitioning_.coordiantesToRank(self_i,self_j+1);  
        int rank_B=partitioning_.coordiantesToRank(self_i,self_j-1);  
        for (int i = 0; i < i_max; i++)
        {
            Buffer_send_B[i]=pressure(i,1);
            Buffer_send_T[i]=pressure(i,j_max-2);
        }
        if (partitioning_.ownRankNo()%2==0) 
        {
          MPI_Send(Buffer_send_T.data(),i_max,MPI_DOUBLE,rank_T,1,MPI_COMM_WORLD);
          MPI_Send(Buffer_send_B.data(),i_max,MPI_DOUBLE,rank_B,1,MPI_COMM_WORLD);           
          MPI_Recv(Buffer_recv_T.data(),j_max,MPI_DOUBLE,rank_T,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
          MPI_Recv(Buffer_recv_B.data(),j_max,MPI_DOUBLE,rank_B,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
       
        } else
        {
          MPI_Recv(Buffer_recv_T.data(),j_max,MPI_DOUBLE,rank_T,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
          MPI_Recv(Buffer_recv_B.data(),j_max,MPI_DOUBLE,rank_B,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
          MPI_Send(Buffer_send_T.data(),i_max,MPI_DOUBLE,rank_T,1,MPI_COMM_WORLD);
          MPI_Send(Buffer_send_B.data(),i_max,MPI_DOUBLE,rank_B,1,MPI_COMM_WORLD);                          
        }
        for (int i = 0; i < i_max; i++)
         {
            pressure(i,j_max-1)=Buffer_recv_T[i];
            pressure(i,0)=Buffer_recv_B[i];
         }        
        }
}  



void Discretization::setBorderVelocityParalell(std::array<double,2> top,std::array<double,2> left,std::array<double,2> right,std::array<double,2> bottom)
{  
   int i_u_max = velocity_X.size()[0], j_u_max = velocity_X.size()[1];
  
   for (int j = 0; j < j_u_max; j++)
   {
        velocity_X(uIBegin(),j)=left[0];
        velocity_X(i_u_max-1,j)=right[0];
   }
   for (int i = 1; i < i_u_max-1; i++) // i starts at 1 and goes to i_u_max-1 so that the wall is the BC in corners
   {
       velocity_X(i,uJBegin())=2*bottom[0]-velocity_X(i,vJBegin()+1);
       velocity_X(i,j_u_max-1)=2*top[0]-velocity_X(i,j_u_max-2);
   }

   // set v velocity
   for (int j = vJBegin(); j < vJEnd(); j++)
   {
       velocity_Y(vIBegin(),j)=2*left[1]-velocity_Y(vIBegin()+1,j); 
       velocity_Y(vIEnd()-1,j)=2*right[1]-velocity_Y(vIEnd()-2,j);

   }
   for (int i = vIBegin()+1; i < vIEnd()-1; i++) // i starts at 1 and goes to i_u_max-1 so that the wall is the BC in corners
   {
       velocity_Y(i,vJBegin())=bottom[1]; 
       velocity_Y(i,vJEnd()-1)=top[1]; 
   }      
}





void Discretization::updateBoundaryFGParalell()
{
    for (int j = uJBegin(); j < uJEnd(); j++)
   {
        F(0,j)=u(0,j);
        F(uIEnd()-1,j)=u(uIEnd()-1,j);
   }

   for (int i = vIBegin(); i < vIEnd(); i++)
   {
        G(i,vJBegin())=v(i,vJBegin());
        G(i,vJEnd()-1)=v(i,vJEnd()-1);

   }
}








