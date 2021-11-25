#pragma once

//#include "computation/computation.h"
#include <cstdlib>
#include <array>

#include "settings.h"



class Partitioning 
{
protected:
        std::array<int,2> offset();
public:
    Partitioning(Settings settings);
    
    int ownRankNo();
    int ownPartitionContainsBottomBoundary();
    int ownPartitionContainsLeftBoundary();
    int ownPartitionContainsRightBoundary();
    int ownPartitionContainsTopBoundary();
    std::array<int,2> nodeOffset();     
};
