#pragma once

#include <cstddef>

#include <astray/parallel/thrust.hpp>

namespace ast
{
template <typename type, bool fortran_order = false>
__device__ constexpr type        unravel_index    (std::size_t index      , const type& dimensions)
{
  type multi_index {};

  if constexpr (fortran_order)
    for (std::size_t i = 0; i < static_cast<std::size_t>(dimensions.size()); ++i)
    {
      multi_index[i] = index % dimensions[i];
      index          = index / dimensions[i];
    }
  else
    for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(dimensions.size() - 1); i >= 0; --i)
    {
      multi_index[i] = index % dimensions[i];
      index          = index / dimensions[i];
    }
  
  return multi_index;
}
template <typename type, bool fortran_order = false>
__device__ constexpr std::size_t ravel_multi_index(const type& multi_index, const type& dimensions)
{
  std::size_t index(0);

  if constexpr (fortran_order)
    for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(dimensions.size() - 1); i >= 0; --i)
      index = index * dimensions[i] + multi_index[i];
  else
    for (std::size_t i = 0; i < static_cast<std::size_t>(dimensions.size()); ++i)
      index = index * dimensions[i] + multi_index[i];
  
  return index;
}
}