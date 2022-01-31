#pragma once

#include <vector>
#include <array>

/** This class represents a 2D array of double values.
 *  Internally they are stored consecutively in memory.
 *  The entries can be accessed by two indices i,j.
 */
class IntArray2D
{
protected:
    std::vector<int> data_;      //< storage array values, in row-major order
    const std::array<int, 2> size_; //< width, height of the domain
public:
    //! constructor
    IntArray2D(std::array<int, 2> size);
    //! destructor
    ~IntArray2D();

    //! get the size
    std::array<int, 2> size() const;

    //! access the value at coordinate (i,j), declared not const, i.e. the value can be changed
    int &operator()(int i, int j);

    //! get the value at coordinate (i,j), declared const, i.e. it is not possible to change the value
    int operator()(int i, int j) const;

    // function returns absolute valkue at position i j
    int abs(int i, int j) const;
};
