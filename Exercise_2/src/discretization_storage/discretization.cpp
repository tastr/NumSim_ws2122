#include "discretization.h"
#include "staggeredgrid.h"
#include "settings.h"
#include "fieldvariable.h"


//constructor
Discretization::Discretization(Settings settings, Partitioning partitioning):StaggeredGrid(settings) //maybe replace settings nCell with partitioning nCells.
,F({settings.nCells[0]+3, settings.nCells[1]+3}, {0,0.5}, {settings.physicalSize[0] / (1.0*settings.nCells[0]), settings.physicalSize[1] / (1.0*settings.nCells[1])}) // is actually smaller than this, but makes handling easier
,G({settings.nCells[0]+3, settings.nCells[1]+3}, {0.5,0}, {settings.physicalSize[0] / (1.0*settings.nCells[0]), settings.physicalSize[1] / (1.0*settings.nCells[1])}) // is actually smaller than this, but makes handling easier
,rhs_({settings.nCells[0]+3, settings.nCells[1]+3}, {0.5,0.5}, {settings.physicalSize[0] / (1.0*settings.nCells[0]), settings.physicalSize[1] / (1.0*settings.nCells[1])}) // is actually smaller than this, but makes handling easier
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
            for (int j = pJBegin(); j < pJEnd(); j++)
            {
            pressure(0,j)=pressure(1,j);
            pressure(i_max-1,j)=pressure(i_max-2,j);
            }
    }else if (partitioning_.ownPartitionContainsLeftBoundary() && !partitioning_.ownPartitionContainsRightBoundary())
        {  std::vector<double> Buffer_send(j_max,0);
        std::vector<double> Buffer_recv(j_max,0);
        int rank=partitioning_.coordiantesToRank(self_i+1,self_j);  
        for (int j = pJBegin(); j < pJEnd(); j++)
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
        for (int j = pJBegin(); j < pJEnd(); j++)
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
          
          for (int j = pJBegin(); j < pJEnd(); j++)
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
        for (int j = pJBegin(); j < pJEnd(); j++)
        {
            Buffer_send_l[j]=pressure(1,j);
            Buffer_send_r[j]=pressure(i_max-2,j);
        }
        if (partitioning_.ownRankNo()%2==0) 
        {
          MPI_Send(Buffer_send_l.data(),j_max,MPI_DOUBLE,rank_l,1,MPI_COMM_WORLD);
          MPI_Send(Buffer_send_r.data(),j_max,MPI_DOUBLE,rank_r,1,MPI_COMM_WORLD);           
          MPI_Recv(Buffer_recv_r.data(),j_max,MPI_DOUBLE,rank_r,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
          MPI_Recv(Buffer_recv_l.data(),j_max,MPI_DOUBLE,rank_l,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
       
        } else
        {
          MPI_Recv(Buffer_recv_r.data(),j_max,MPI_DOUBLE,rank_r,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
          MPI_Recv(Buffer_recv_l.data(),j_max,MPI_DOUBLE,rank_l,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
          MPI_Send(Buffer_send_l.data(),j_max,MPI_DOUBLE,rank_l,1,MPI_COMM_WORLD);
          MPI_Send(Buffer_send_r.data(),j_max,MPI_DOUBLE,rank_r,1,MPI_COMM_WORLD);               
        }
        for (int j = pJBegin(); j < pJEnd(); j++)
         {
            pressure(0,j)=Buffer_recv_l[j];
            pressure(i_max-1,j)=Buffer_recv_r[j];
         }        
        }
////////////////////////////////////////////////////////////////////////////////
   if (partitioning_.ownPartitionContainsTopBoundary() && partitioning_.ownPartitionContainsBottomBoundary())
       for (int i = pIBegin(); i < pIEnd(); i++)
        {
            pressure(i,0)=pressure(i,1);
            pressure(i,j_max-1)=pressure(i,j_max-2);
        }
    if (!partitioning_.ownPartitionContainsTopBoundary() && partitioning_.ownPartitionContainsBottomBoundary())
       { 
        std::vector<double> Buffer_send(i_max,0);
        std::vector<double> Buffer_recv(i_max,0);
        int rank=partitioning_.coordiantesToRank(self_i,self_j+1);    
        for (int i = pIBegin(); i < pIEnd(); i++)
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
        for (int i = pIBegin(); i < pIEnd(); i++)
            {
                pressure(i,j_max-1)=Buffer_recv[i];
            }
        }    
        if (partitioning_.ownPartitionContainsTopBoundary() && !partitioning_.ownPartitionContainsBottomBoundary())
        {
        std::vector<double> Buffer_send(i_max,0);
        std::vector<double> Buffer_recv(i_max,0);
        int rank=partitioning_.coordiantesToRank(self_i,self_j-1);    
        for (int i = pIBegin(); i < pIEnd(); i++)
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
        for (int i = pIBegin(); i < pIEnd(); i++)
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
        for (int i = pIBegin(); i < pIEnd(); i++)
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
        for (int i = pIBegin(); i < pIEnd(); i++)
         {
            pressure(i,j_max-1)=Buffer_recv_T[i];
            pressure(i,0)=Buffer_recv_B[i];
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
   int lu=(j_u_max-2)-(uJBegin()+1);
   int lv=(j_v_max-2)-(vJBegin()+1);
   
   if (partitioning_.ownPartitionContainsLeftBoundary() && partitioning_.ownPartitionContainsRightBoundary())
   for (int j = uJBegin()+1; j < j_u_max-2; j++)
   {
        velocity_X(i_u_begin+1,j)=left[0];
        velocity_X(i_u_max-2,j)=right[0];
   }
    for (int j = vJBegin(); j < j_v_max-2; j++)
   {
       velocity_Y(i_v_begin+2,j)=2*left[1]-velocity_Y(i_v_begin+1,j); 
       velocity_Y(i_v_max-1,j)=2*right[1]-velocity_Y(i_v_max-2,j);
   }

   if (partitioning_.ownPartitionContainsLeftBoundary() && !partitioning_.ownPartitionContainsRightBoundary())
   {    std::vector<double> Buffer_send_u(2*lu,0);
        std::vector<double> Buffer_send_v(lv,0);
        std::vector<double> Buffer_recv_u(lu,0);
        std::vector<double> Buffer_recv_v(lv,0);
        int rank=partitioning_.coordiantesToRank(self_i+1,self_j);  

        for (int j = uJBegin()+1; j < j_u_max-2; j++)
           {
            velocity_X(i_u_begin,j)=left[0];
            Buffer_send_u[j-(j_u_begin+1)]=velocity_X(i_u_max-2,j);
            Buffer_send_u[j_u_max-2+j-(j_u_begin+1)]=velocity_X(i_u_max-3,j);
           }
        for (int j = vJBegin()+1; j < j_v_max-2; j++)
           {
            velocity_Y(i_v_begin+2,j)=2*left[1]-velocity_Y(i_v_begin+1,j); 
            Buffer_send_v[j-(vJBegin()+1)]= velocity_Y(i_v_max-2,j);
           }
        if (partitioning_.ownRankNo()%2==0) 
        {
          MPI_Send(Buffer_send_u.data(),2*lu,MPI_DOUBLE,rank,1,MPI_COMM_WORLD);
          MPI_Send(Buffer_send_v.data(),lv,MPI_DOUBLE,rank,1,MPI_COMM_WORLD);           
          MPI_Recv(Buffer_recv_u.data(),lu,MPI_DOUBLE,rank,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
          MPI_Recv(Buffer_recv_v.data(),lv,MPI_DOUBLE,rank,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
        }else
        {
          MPI_Recv(Buffer_recv_u.data(),lu,MPI_DOUBLE,rank,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
          MPI_Recv(Buffer_recv_v.data(),lv,MPI_DOUBLE,rank,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
          MPI_Send(Buffer_send_u.data(),2*lu,MPI_DOUBLE,rank,1,MPI_COMM_WORLD);
          MPI_Send(Buffer_send_v.data(),lv,MPI_DOUBLE,rank,1,MPI_COMM_WORLD);            
        }
        for (int j = uJBegin()+1; j < j_u_max-2; j++)
        {
            velocity_X(i_u_max-1,j)=Buffer_recv_u[j-(uJBegin()+1)];
        } 
        for (int j = vJBegin()+1; j < j_v_max-2; j++)
         {
            velocity_Y(i_v_max-1,j)=Buffer_recv_v[j-(vJBegin()+1)];
        }
    }
    if (!partitioning_.ownPartitionContainsLeftBoundary() && partitioning_.ownPartitionContainsRightBoundary())
   {    std::vector<double> Buffer_send_u(lu,0);
        std::vector<double> Buffer_send_v(lv,0);
        std::vector<double> Buffer_recv_u(2*lu,0);
        std::vector<double> Buffer_recv_v(lv,0);
        int rank=partitioning_.coordiantesToRank(self_i-1,self_j);  
        for (int j = uJBegin()+1; j < j_u_max-2; j++)
           {
            velocity_X(i_u_max-2,j)=left[0];
            Buffer_send_u[j-(j_u_begin+1)]=velocity_X(i_u_begin+2,j);
            }
        for (int j = vJBegin()+1; j < j_v_max-2; j++)
           {
            velocity_Y(i_v_max-1,j)=2*left[1]-velocity_Y(i_v_max-2,j); 
            Buffer_send_v[j-(vJBegin()+1)]=velocity_Y(i_v_begin+1,j);
           }
        if (partitioning_.ownRankNo()%2==0) 
        {
          MPI_Send(Buffer_send_u.data(),lu,MPI_DOUBLE,rank,1,MPI_COMM_WORLD);
          MPI_Send(Buffer_send_v.data(),lv,MPI_DOUBLE,rank,1,MPI_COMM_WORLD);           
          MPI_Recv(Buffer_recv_u.data(),2*lu,MPI_DOUBLE,rank,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
          MPI_Recv(Buffer_recv_v.data(),lv,MPI_DOUBLE,rank,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
        }else
        {
          MPI_Recv(Buffer_recv_u.data(),2*lu,MPI_DOUBLE,rank,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
          MPI_Recv(Buffer_recv_v.data(),lv,MPI_DOUBLE,rank,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
          MPI_Send(Buffer_send_u.data(),lu,MPI_DOUBLE,rank,1,MPI_COMM_WORLD);
          MPI_Send(Buffer_send_v.data(),lv,MPI_DOUBLE,rank,1,MPI_COMM_WORLD);            
        }
        for (int j = uJBegin()+1; j < j_u_max-2; j++)
        {
            velocity_X(i_u_begin+1,j)=Buffer_recv_u[j-(uJBegin()+1)];
            velocity_X(i_u_begin,j)=Buffer_recv_u[(j_u_max-2)+j-(uJBegin()+1)];
        } 
        for (int j = vJBegin()+1; j < j_v_max-2; j++)
         {
            velocity_Y(i_v_begin,j)=Buffer_recv_v[j-(vJBegin()+1)];
        }
    }else
    {
        std::vector<double> Buffer_send_u_l(lu,0);
        std::vector<double> Buffer_send_v_l(lv,0);
        std::vector<double> Buffer_recv_u_l(2*lu,0);
        std::vector<double> Buffer_recv_v_l(lv,0);
        std::vector<double> Buffer_send_u_r(2*lu,0);
        std::vector<double> Buffer_send_v_r(lv,0);
        std::vector<double> Buffer_recv_u_r(2*lu,0);
        std::vector<double> Buffer_recv_v_r(lv,0);
        int rank_l=partitioning_.coordiantesToRank(self_i-1,self_j);  
        int rank_r=partitioning_.coordiantesToRank(self_i+1,self_j);  
        for (int j = uJBegin()+1; j < j_u_max-2; j++)
           {
            Buffer_send_u_r[j-(j_u_begin+1)]=velocity_X(i_u_begin+2,j);
            Buffer_send_u_l[j-(j_u_begin+1)]=velocity_X(i_u_max-2,j);
            Buffer_send_u_l[j_u_max-2+j-(j_u_begin+1)]=velocity_X(i_u_max-3,j);
        
            }
        for (int j = vJBegin()+1; j < j_v_max-2; j++)
           {
            Buffer_send_v_r[j-(vJBegin()+1)]=velocity_Y(i_v_begin+1,j);
            Buffer_send_v_l[j-(vJBegin()+1)]= velocity_Y(i_v_max-2,j);
           }
        if (partitioning_.ownRankNo()%2==0) 
        {
          MPI_Send(Buffer_send_u_r.data(),2*lu,MPI_DOUBLE,rank_r,1,MPI_COMM_WORLD);
          MPI_Send(Buffer_send_v_r.data(),lv,MPI_DOUBLE,rank_r,1,MPI_COMM_WORLD);           
          MPI_Send(Buffer_send_u_l.data(),lu,MPI_DOUBLE,rank_l,1,MPI_COMM_WORLD);
          MPI_Send(Buffer_send_v_l.data(),lv,MPI_DOUBLE,rank_l,1,MPI_COMM_WORLD);           
          MPI_Recv(Buffer_recv_u_r.data(),lu,MPI_DOUBLE,rank_r,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
          MPI_Recv(Buffer_recv_v_r.data(),lv,MPI_DOUBLE,rank_r,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
          MPI_Recv(Buffer_recv_u_l.data(),2*lu,MPI_DOUBLE,rank_l,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
          MPI_Recv(Buffer_recv_v_l.data(),lv,MPI_DOUBLE,rank_l,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
        }else
        {
          MPI_Recv(Buffer_recv_u_r.data(),lu,MPI_DOUBLE,rank_r,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
          MPI_Recv(Buffer_recv_v_r.data(),lv,MPI_DOUBLE,rank_r,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
          MPI_Recv(Buffer_recv_u_l.data(),2*lu,MPI_DOUBLE,rank_l,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
          MPI_Recv(Buffer_recv_v_l.data(),lv,MPI_DOUBLE,rank_l,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
          MPI_Send(Buffer_send_u_r.data(),2*lu,MPI_DOUBLE,rank_r,1,MPI_COMM_WORLD);
          MPI_Send(Buffer_send_v_r.data(),lv,MPI_DOUBLE,rank_r,1,MPI_COMM_WORLD);           
          MPI_Send(Buffer_send_u_l.data(),lu,MPI_DOUBLE,rank_l,1,MPI_COMM_WORLD);
          MPI_Send(Buffer_send_v_l.data(),lv,MPI_DOUBLE,rank_l,1,MPI_COMM_WORLD);           
        }   
            for (int j = uJBegin()+1; j < j_u_max-2; j++)
        {
            velocity_X(i_u_begin+1,j)=Buffer_recv_u_l[j-(uJBegin()+1)];
            velocity_X(i_u_begin,j)=Buffer_recv_u_l[(j_u_max-2)+j-(uJBegin()+1)];
            velocity_X(i_u_max-1,j)=Buffer_recv_u_r[j-(uJBegin()+1)];
        
        } 
        for (int j = vJBegin()+1; j < j_v_max-2; j++)
         {
            velocity_Y(i_v_begin,j)=Buffer_recv_v_l[j-(vJBegin()+1)];
            velocity_Y(i_v_max-1,j)=Buffer_recv_v_r[j-(vJBegin()+1)];
        
        }
    }
    // top bottom
    lu=(uIEnd()-2)-(uIBegin()+1);
    lv=(vIEnd()-2)-(vIBegin()+1);
    if (partitioning_.ownPartitionContainsTopBoundary() && partitioning_.ownPartitionContainsBottomBoundary())
    {
        for (int i = uIBegin()+1; i < uIEnd()-2; i++) // i starts at 1 and goes to i_u_max-1 so that the wall is the BC in corners
         {
         velocity_X(i,j_u_begin)=2*bottom[0]-velocity_X(i,j_u_begin+1);
         velocity_X(i,j_u_max-2)=2*top[0]-velocity_X(i,j_u_max-3);
         }
       for (int i = vIBegin()+1; i < vIEnd()-2; i++) // i starts at 1 and goes to i_u_max-1 so that the wall is the BC in corners
        {
         velocity_Y(i,j_v_begin)=bottom[1]; 
         velocity_Y(i,j_v_max-2)=top[1]; 
         } 
    }else  if (!partitioning_.ownPartitionContainsTopBoundary() && partitioning_.ownPartitionContainsBottomBoundary())
   {    
        std::vector<double> Buffer_send_u_T(lu,0);
        std::vector<double> Buffer_send_v_T(lv,0);
        std::vector<double> Buffer_recv_u_T(lu,0);
        std::vector<double> Buffer_recv_v_T(2*lv,0);
        int rank_T=partitioning_.coordiantesToRank(self_i,self_j+1);  
        
        for (int i = uIBegin()+1; i < uIEnd()-2; i++) // i starts at 1 and goes to i_u_max-1 so that the wall is the BC in corners
         {
         velocity_X(i,j_u_begin)=2*bottom[0]-velocity_X(i,j_u_begin+1);
         Buffer_send_u_T[i-(uIBegin()+1)]=velocity_X(i,j_u_max-3); //hier hab ich nen anderen index, -2 waere doch der wert im ghostlayer
         }
        for (int i = vIBegin()+1; i < vIEnd()-2; i++) // i starts at 1 and goes to i_u_max-1 so that the wall is the BC in corners
         {
         velocity_Y(i,j_v_begin)=bottom[1]; 
         Buffer_send_v_T[i-(vIBegin()+1)]=velocity_Y(i,j_v_max-3);
         }
        if (partitioning_.ownRankNo()%2==0) 
        {
          MPI_Send(Buffer_send_u_T.data(),lu,MPI_DOUBLE,rank_T,1,MPI_COMM_WORLD);
          MPI_Send(Buffer_send_v_T.data(),lv,MPI_DOUBLE,rank_T,1,MPI_COMM_WORLD);           
          MPI_Recv(Buffer_recv_u_T.data(),lu,MPI_DOUBLE,rank_T,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
          MPI_Recv(Buffer_recv_v_T.data(),2*lv,MPI_DOUBLE,rank_T,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
        }else
        {
          MPI_Recv(Buffer_recv_u_T.data(),lu,MPI_DOUBLE,rank_T,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
          MPI_Recv(Buffer_recv_v_T.data(),2*lv,MPI_DOUBLE,rank_T,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
          MPI_Send(Buffer_send_u_T.data(),lu,MPI_DOUBLE,rank_T,1,MPI_COMM_WORLD);
          MPI_Send(Buffer_send_v_T.data(),lv,MPI_DOUBLE,rank_T,1,MPI_COMM_WORLD);            
        }
        for (int i = uIBegin()+1; i < uIEnd()-2; i++) // i starts at 1 and goes to i_u_max-1 so that the wall is the BC in corners
         {
         velocity_X(i,j_u_max-2)=Buffer_recv_u_T[i-(uIBegin()+1)];
         }
        for (int i = vIBegin()+1; i < vIEnd()-2; i++) // i starts at 1 and goes to i_u_max-1 so that the wall is the BC in corners
         {
         velocity_Y(i,j_v_max-1)=Buffer_recv_v_T[i-(vIBegin()+1)];
         velocity_Y(i,j_v_max-2)=Buffer_recv_v_T[(vIEnd()-2)+i-(vIBegin()+1)]; 
         }
    }else  if (partitioning_.ownPartitionContainsTopBoundary() && !partitioning_.ownPartitionContainsBottomBoundary())
     { 
        std::vector<double> Buffer_send_u_B(lu,0);
        std::vector<double> Buffer_send_v_B(2*lv,0);
        std::vector<double> Buffer_recv_u_B(lu,0);
        std::vector<double> Buffer_recv_v_B(lv,0);
        int rank_B=partitioning_.coordiantesToRank(self_i,self_j-1);  
        for (int i = uIBegin()+1; i < uIEnd()-2; i++) // i starts at 1 and goes to i_u_max-1 so that the wall is the BC in corners
         {
         velocity_X(i,j_u_max-2)=2*top[0]-velocity_X(i,j_u_max-3);
         Buffer_send_u_B[i-(uIBegin()+1)]= velocity_X(i,j_u_begin+1);
         }
        for (int i = vIBegin()+1; i < vIEnd()-2; i++) // i starts at 1 and goes to i_u_max-1 so that the wall is the BC in corners
        {
         velocity_Y(i,j_v_max-2)=top[1]; 
         Buffer_send_v_B[i-(vIBegin()+1)]=velocity_Y(i,j_v_begin+3);
         Buffer_send_v_B[(vIEnd()-2)+i-(vIBegin()+1)]=velocity_Y(i,j_v_begin+2);
        
        }
        if (partitioning_.ownRankNo()%2==0) 
        {
          MPI_Send(Buffer_send_u_B.data(),lu,MPI_DOUBLE,rank_B,1,MPI_COMM_WORLD);
          MPI_Send(Buffer_send_v_B.data(),2*lv,MPI_DOUBLE,rank_B,1,MPI_COMM_WORLD);           
          MPI_Recv(Buffer_recv_u_B.data(),lu,MPI_DOUBLE,rank_B,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
          MPI_Recv(Buffer_recv_v_B.data(),lv,MPI_DOUBLE,rank_B,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
        }else
        {
          MPI_Recv(Buffer_recv_u_B.data(),lu,MPI_DOUBLE,rank_B,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
          MPI_Recv(Buffer_recv_v_B.data(),lv,MPI_DOUBLE,rank_B,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
          MPI_Send(Buffer_send_u_B.data(),lu,MPI_DOUBLE,rank_B,1,MPI_COMM_WORLD);
          MPI_Send(Buffer_send_v_B.data(),2*lv,MPI_DOUBLE,rank_B,1,MPI_COMM_WORLD);            
        }
         for (int i = uIBegin()+1; i < uIEnd()-2; i++) // i starts at 1 and goes to i_u_max-1 so that the wall is the BC in corners
         {
         velocity_X(i,j_u_begin)=Buffer_recv_u_B[i-(uIBegin()+1)];
         }
        for (int i = vIBegin()+1; i < vIEnd()-2; i++) // i starts at 1 and goes to i_u_max-1 so that the wall is the BC in corners
         {
         velocity_Y(i,j_v_begin)=Buffer_recv_v_B[i-(vIBegin()+1)];
         }
    }else
    {   std::vector<double> Buffer_send_u_T(lu,0);
        std::vector<double> Buffer_send_v_T(lv,0);
        std::vector<double> Buffer_recv_u_T(lu,0);
        std::vector<double> Buffer_recv_v_T(2*lv,0);
        std::vector<double> Buffer_send_u_B(lu,0);
        std::vector<double> Buffer_send_v_B(2*lv,0);
        std::vector<double> Buffer_recv_u_B(lu,0);
        std::vector<double> Buffer_recv_v_B(lv,0);
        int rank_B=partitioning_.coordiantesToRank(self_i,self_j-1);  
        int rank_T=partitioning_.coordiantesToRank(self_i,self_j+1);  
        
        for (int i = uIBegin()+1; i < uIEnd()-2; i++) // i starts at 1 and goes to i_u_max-1 so that the wall is the BC in corners
         {
          Buffer_send_u_B[i-(uIBegin()+1)]= velocity_X(i,j_u_begin+1);         
          Buffer_send_u_T[i-(uIBegin()+1)]=velocity_X(i,j_u_max-3); //hier hab ich nen anderen index, -2 waere doch der wert im ghostlayer
         }
        for (int i = vIBegin()+1; i < vIEnd()-2; i++) // i starts at 1 and goes to i_u_max-1 so that the wall is the BC in corners
         {
         velocity_Y(i,j_v_begin)=bottom[1]; 
         Buffer_send_v_T[i-(vIBegin()+1)]=velocity_Y(i,j_v_max-3);
         Buffer_send_v_B[i-(vIBegin()+1)]=velocity_Y(i,j_v_begin+3);
         Buffer_send_v_B[(vIEnd()-2)+i-(vIBegin()+1)]=velocity_Y(i,j_v_begin+2);
        
         }
        if (partitioning_.ownRankNo()%2==0) 
        {
          MPI_Send(Buffer_send_u_T.data(),lu,MPI_DOUBLE,rank_T,1,MPI_COMM_WORLD);
          MPI_Send(Buffer_send_v_T.data(),lv,MPI_DOUBLE,rank_T,1,MPI_COMM_WORLD);           
          MPI_Send(Buffer_send_u_B.data(),lu,MPI_DOUBLE,rank_B,1,MPI_COMM_WORLD);
          MPI_Send(Buffer_send_v_B.data(),2*lv,MPI_DOUBLE,rank_B,1,MPI_COMM_WORLD);           
          MPI_Recv(Buffer_recv_u_T.data(),lu,MPI_DOUBLE,rank_T,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
          MPI_Recv(Buffer_recv_v_T.data(),2*lv,MPI_DOUBLE,rank_T,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
          MPI_Recv(Buffer_recv_u_B.data(),lu,MPI_DOUBLE,rank_B,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
          MPI_Recv(Buffer_recv_v_B.data(),lv,MPI_DOUBLE,rank_B,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
        }else
        {
          MPI_Recv(Buffer_recv_u_T.data(),lu,MPI_DOUBLE,rank_T,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
          MPI_Recv(Buffer_recv_v_T.data(),2*lv,MPI_DOUBLE,rank_T,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);    
          MPI_Recv(Buffer_recv_u_B.data(),lu,MPI_DOUBLE,rank_B,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE); 
          MPI_Recv(Buffer_recv_v_B.data(),lv,MPI_DOUBLE,rank_B,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);           
          MPI_Send(Buffer_send_u_T.data(),lu,MPI_DOUBLE,rank_T,1,MPI_COMM_WORLD);
          MPI_Send(Buffer_send_v_T.data(),lv,MPI_DOUBLE,rank_T,1,MPI_COMM_WORLD);            
          MPI_Send(Buffer_send_u_B.data(),lu,MPI_DOUBLE,rank_B,1,MPI_COMM_WORLD);
          MPI_Send(Buffer_send_v_B.data(),2*lv,MPI_DOUBLE,rank_B,1,MPI_COMM_WORLD);           
          
        }
        for (int i = uIBegin()+1; i < uIEnd()-2; i++) // i starts at 1 and goes to i_u_max-1 so that the wall is the BC in corners
         {
         velocity_X(i,j_u_max-2)=Buffer_recv_u_T[i-(uIBegin()+1)];
         velocity_X(i,j_u_begin)=Buffer_recv_u_B[i-(uIBegin()+1)];
         
         }
        for (int i = vIBegin()+1; i < vIEnd()-2; i++) // i starts at 1 and goes to i_u_max-1 so that the wall is the BC in corners
         {
         velocity_Y(i,j_v_max-1)=Buffer_recv_v_T[i-(vIBegin()+1)];
         velocity_Y(i,j_v_max-2)=Buffer_recv_v_T[(vIEnd()-2)+i-(vIBegin()+1)]; 
         velocity_Y(i,j_v_begin)=Buffer_recv_v_B[i-(vIBegin()+1)];
         
        }
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








