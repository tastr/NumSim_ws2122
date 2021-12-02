#pragma once

//#include "computation/computation.h"
#include <cstdlib>
#include <array>
#include "settings.h"
#include <mpi.h>



class Partitioning 
{
protected:
        // std::array<int,2> offset(); // Scheint nicht implementiert zu sein
        int n; // Number of Subdomains horizontally
        int m; // Number of Subdomains verticallz
        
    int ownPartitionContainsBottomBoundaryValue;
    int ownPartitionContainsLeftBoundaryValue;
    int ownPartitionContainsRightBoundaryValue;
    int ownPartitionContainsTopBoundaryValue;
    int ownRankNoValue;
    int firstValue;
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
    void setFirst();
    void setSubgrid();

public:
    Partitioning(Settings settings);
    
    int ownRankNo() const;
    int ownPartitionContainsBottomBoundary() const;
    int ownPartitionContainsLeftBoundary() const;
    int ownPartitionContainsRightBoundary() const;
    int ownPartitionContainsTopBoundary() const;
    std::array<int,2> nodeOffset() const; // returns how many nodes are before this node in each direction
    std::array<int,2> cellsNodeOffset() const; // returns how many cells this node is offset
    std::array<int,2> nCells() const;
    std::array<int,2> nCellsGlobal() const; 
    int coordiantesToRank(int i, int j) const;     
    int first() const;
    std::array<int,2> getSubGridSize() const; 
};
