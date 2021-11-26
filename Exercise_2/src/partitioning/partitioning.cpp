#include "partitioning.h"

Partitioning::Partitioning(Settings settings):
nCellsGlobal_(settings.nCells)
{
    n=3;  // predefined for Tests, use 3x3 since it contains all possible bordercases 
    m=3;  // in the project these numbers need to be calculated in a function that defines the domainsplitting.
    setOwnRankNo();//could be given to the partitioning
    setNodeOffset();
    setOwnPartitionContainsBottomBoundary();
    setOwnPartitionContainsLeftBoundary();
    setOwnPartitionContainsRightBoundary();
    setOwnPartitionContainsTopBoundary();  
    setNCells();
    std::cout<< n << m <<std::endl; 
    std::cout<< ownRankNoValue <<std::endl;  
    std::cout<< nodeOffsetValue[0] << nodeOffsetValue[1] << std::endl;  
    std::cout<<  nCells[0] <<  nCells[1] << std::endl;
    }




    int Partitioning::ownRankNo()
    {
        return ownRankNoValue;

    }
    int Partitioning::ownPartitionContainsBottomBoundary()
    {
         return ownPartitionContainsBottomBoundaryValue;
    }
    int Partitioning::ownPartitionContainsLeftBoundary()
    {
         return ownPartitionContainsLeftBoundaryValue;
    }
    int Partitioning::ownPartitionContainsRightBoundary()
    {
         return ownPartitionContainsRightBoundaryValue;
    }
    int Partitioning::ownPartitionContainsTopBoundary()
    {
         return ownPartitionContainsTopBoundaryValue;
    }
    std::array<int,2> Partitioning::nodeOffset()
    {
       return nodeOffsetValue;
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
    //return MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    ownRankNoValue=8;
    }

    void Partitioning::setNodeOffset()
    {int i=0;
     int j=0;
     
    while(ownRankNoValue-(j+1)*n>0)
        {        
        j++;
        } 
    i=ownRankNoValue-j*n;
    nodeOffsetValue={i,j};  
    //enumeration from 0 to n-1/m-1
    }

    void Partitioning::setNCells()
    {
        nCells[0]=nCellsGlobal_[0]/n;
        nCells[1]=nCellsGlobal_[1]/m;
        if (nodeOffsetValue[0]==n-1)
        {
            nCells[0]=nCellsGlobal_[0]-(nodeOffsetValue[0])*nCells[0];
        }
        if (nodeOffsetValue[1]==m-1)
        {
            nCells[1]=nCellsGlobal_[1]-(nodeOffsetValue[1])*nCells[1];
        }
         
   
    }
















