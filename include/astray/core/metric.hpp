#pragma once

#include <astray/core/termination_reason.hpp>
#include <astray/math/coordinate_system.hpp>
#include <astray/math/linear_algebra.hpp>
#include <astray/parallel/thrust.hpp>

namespace ast
{
template <
  coordinate_system_type system                   ,
  typename               scalar_type              , 
  typename               vector_type              = vector4  <scalar_type>, 
  typename               christoffel_symbols_type = tensor444<scalar_type>>
class metric
{
public:
  virtual ~metric() = default;

  static constexpr coordinate_system_type     coordinate_system          ()
  {
    return system;
  }
  
  __device__ virtual scalar_type              coordinate_system_parameter() const
  {
    return scalar_type(0);
  }
  __device__ virtual termination_reason       check_termination          (const vector_type& position, const vector_type& direction) const
  {
    return termination_reason::none;
  }
  
  __device__ virtual christoffel_symbols_type christoffel_symbols        (const vector_type& position) const = 0;
};
}