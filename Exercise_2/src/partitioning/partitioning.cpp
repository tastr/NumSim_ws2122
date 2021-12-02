#include "partitioning.h"

Partitioning::Partitioning(Settings settings):
nCellsGlobal_(settings.nCells)
{   // n=2 m=3 produziert error
    //n=1;  // predefined for Tests, use 3x3 since it contains all possible bordercases 
    //m=2;  // in the project these numbers need to be calculated in a function that defines the domainsplitting.
    setSubgrid();
    setOwnRankNo();//could be given to the partitioning
    setNodeOffset();
    setOwnPartitionContainsBottomBoundary();
    setOwnPartitionContainsLeftBoundary();
    setOwnPartitionContainsRightBoundary();
    setOwnPartitionContainsTopBoundary();  
    setNCells();
    setFirst();

    //std::cout<< n << m <<std::endl; 
    //std::cout<< ownRankNoValue <<std::endl;  
    //std::cout<<"Node offset "<< nodeOffsetValue[0] << nodeOffsetValue[1] << std::endl;  
    //std::cout<<  nCells_[0] <<  nCells_[1] << std::endl;
    //std::cout<< ownRankNo() << " Left "<< ownPartitionContainsLeftBoundary() <<std::endl;
    //std::cout<< ownRankNo() << " Right "<< ownPartitionContainsRightBoundary() <<std::endl;
    //std::cout<< ownRankNo() << " TOP "<< ownPartitionContainsTopBoundary() <<std::endl;
    //std::cout<< ownRankNo() << " Bottom "<< ownPartitionContainsBottomBoundary() <<std::endl;
   
}




    int Partitioning::ownRankNo() const
    {
        return ownRankNoValue;

    }
    int Partitioning::ownPartitionContainsBottomBoundary() const
    {
         return ownPartitionContainsBottomBoundaryValue;
    }
    int Partitioning::ownPartitionContainsLeftBoundary() const
    {
         return ownPartitionContainsLeftBoundaryValue;
    }
    int Partitioning::ownPartitionContainsRightBoundary() const
    {
         return ownPartitionContainsRightBoundaryValue;
    }
    int Partitioning::ownPartitionContainsTopBoundary() const
    {
         return ownPartitionContainsTopBoundaryValue;
    }
    std::array<int,2> Partitioning::nodeOffset() const
    {
       return nodeOffsetValue;
    }

   std::array<int,2> Partitioning::nCells() const
    {
     return nCells_;
    } 
    std::array<int,2> Partitioning::nCellsGlobal() const
    {
     return nCellsGlobal_;
    } 
    
    int Partitioning::coordiantesToRank(int i, int j) const
    {
        return  j*n + i;
    }


















    // Initialization functions
    void Partitioning::setOwnPartitionContainsBottomBoundary()
    {
        int result=0;
        if (nodeOffsetValue[1]==0)
        {
            result=1;
        }
        ownPartitionContainsBottomBoundaryValue=result;
    }
    void Partitioning::setOwnPartitionContainsLeftBoundary()
    {
        int result=0;
        if (nodeOffsetValue[0]==0)
        {
            result=1;
        }
        ownPartitionContainsLeftBoundaryValue=result;
    }
    void Partitioning::setOwnPartitionContainsRightBoundary()
    {
        int result=0;
        if (nodeOffsetValue[0]==n-1)
        {
            result=1;
        }
        ownPartitionContainsRightBoundaryValue=result;
    }
    void Partitioning::setOwnPartitionContainsTopBoundary()
    {
        int result=0;
        if (nodeOffsetValue[1]==m-1)
        {
            result=1;
        }
        ownPartitionContainsTopBoundaryValue=result;
    }
    void Partitioning::setOwnRankNo()  
    {
     MPI_Comm_rank(MPI_COMM_WORLD, &ownRankNoValue);
    //ownRankNoValue=7;
    }

    void Partitioning::setNodeOffset()
    {int i=0;
     int j=0;
    //if (n>1)
   // {
    while((ownRankNoValue+1)-(j+1)*n>0)
        {        
        j++;
       } 
    i=ownRankNoValue-j*n;
    //enumeration from 0 to n-1/m-1
   // }else
   // {
   //     i=0;
    //    j=ownRankNoValue;
   // }
    nodeOffsetValue={i,j};  
    }

    void Partitioning::setNCells()
    {
        nCells_[0]=nCellsGlobal_[0]/n;
        nCells_[1]=nCellsGlobal_[1]/m;
        if (nodeOffsetValue[0]==n-1)
        {
            nCells_[0]=nCellsGlobal_[0]-(nodeOffsetValue[0])*nCells_[0];
        }
        if (nodeOffsetValue[1]==m-1)
        {
            nCells_[1]=nCellsGlobal_[1]-(nodeOffsetValue[1])*nCells_[1];
        }
         
   
    }


    int Partitioning::first() const
    {
        return firstValue;
    }
    

    void Partitioning::setFirst()
    {
        firstValue=0;
        if ((nodeOffset()[0]%2==0 && nodeOffset()[1]%2==0) || (nodeOffset()[0]%2==1 && nodeOffset()[1]%2==1) )
        {
            firstValue=1;
        }
    }
 




 void Partitioning::setSubgrid()
{
 int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  n=world_size;
  m=1; 

 int Xcells=nCellsGlobal_[0]/n;
 int Ycells=nCellsGlobal_[1]/m;

int newXcells;
int newYcells;
  
  for (int i = 1; i < world_size+1; i++)
  {
    for (int j = 1; j < world_size+1 ; j++)
    {
     if (i*j == world_size)
    {
        newXcells=nCellsGlobal_[0]/i;
        newYcells=nCellsGlobal_[1]/j;
        if ((newXcells+newYcells) <( Ycells+Xcells))
        {
           n=i;
           m=j;
           Xcells=nCellsGlobal_[0]/n;
           Ycells=nCellsGlobal_[1]/m;
        }
    }

    }
    
  }
  
}


std::array<int,2> Partitioning::getSubGridSize() const
    {   
     return {n,m};
    } 
   







