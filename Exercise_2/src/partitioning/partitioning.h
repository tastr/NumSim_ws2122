#pragma once

//#include "computation/computation.h"
#include <cstdlib>
#include <array>
#include "settings.h"
#include <mpi.h>



class Partitioning 
{
protected:
        std::array<int,2> offset();
        int n; // Number of Subdomains horizontally
        int m; // Number of Subdomains verticallz
        
    int ownPartitionContainsBottomBoundaryValue;
    int ownPartitionContainsLeftBoundaryValue;
    int ownPartitionContainsRightBoundaryValue;
    int ownPartitionContainsTopBoundaryValue;
    int ownRankNoValue;
    std::array<int,2> nodeOffsetValue;
    std::array<int,2> nCellsGlobal_;
    std::array<int,2> nCells_;

    // Initialization functions
    void setOwnPartitionContainsBottomBoundary();
    void setOwnPartitionContainsLeftBoundary();
    void setOwnPartitionContainsRightBoundary();
    void setOwnPartitionContainsTopBoundary();
    void setOwnRankNo(); 
    void setNodeOffset(); 
    void setNCells();

public:
    Partitioning(Settings settings);
    
    int ownRankNo();
    int ownPartitionContainsBottomBoundary();
    int ownPartitionContainsLeftBoundary();
    int ownPartitionContainsRightBoundary();
    int ownPartitionContainsTopBoundary();
    std::array<int,2> nodeOffset(); 
    std::array<int,2> nCells(); 
    std::array<int,2> nCellsGlobal(); 
    int coordiantesToRank(int i, int j);     
  
};
