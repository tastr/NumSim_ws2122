#include "IntArray2d.h"
#include <cassert>

IntArray2D::IntArray2D(std::array<int, 2> size) : size_(size)
{
  // allocate data, initialize to 0
  data_.resize(size_[0] * size_[1], 0.0);
}

IntArray2D::~IntArray2D()
{
}

//! get the size
std::array<int, 2> IntArray2D::size() const
{
  return size_;
}

int &IntArray2D::operator()(int i, int j)
{
  const int index = j * size_[0] + i;

  // assert that indices are in range
  assert(0 <= i && i < size_[0]);
  assert(0 <= j && j < size_[1]);
  assert(j * size_[0] + i < (int)data_.size());

  return data_[index];
}

int IntArray2D::operator()(int i, int j) const
{
  const int index = j * size_[0] + i;

  // assert that indices are in range
  assert(0 <= i && i < size_[0]);
  assert(0 <= j && j < size_[1]);
  assert(j * size_[0] + i < (int)data_.size());

  return data_[index];
}

int IntArray2D::abs(int i, int j) const
{
  const int index = j * size_[0] + i;

  // assert that indices are in range
  assert(0 <= i && i < size_[0]);
  assert(0 <= j && j < size_[1]);
  assert(j * size_[0] + i < (int)data_.size());
  int value = data_[index];

  if (value >= 0)
  {
    return value;
  }
  else
  {
    return -value;
  }
}