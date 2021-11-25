#include "partitioning.h"

Partitioning::Partitioning(Settings settings)
{
    
}


    int Partitioning::ownRankNo()
    {
        return 0;
    }
    int Partitioning::ownPartitionContainsBottomBoundary()
    {
         return 0;
    }
    int Partitioning::ownPartitionContainsLeftBoundary()
    {
         return 0;
    }
    int Partitioning::ownPartitionContainsRightBoundary()
    {
         return 0;
    }
    int Partitioning::ownPartitionContainsTopBoundary()
    {
         return 0;
    }
    std::array<int,2> nodeOffset()
    {
        return {0,0};
    }





