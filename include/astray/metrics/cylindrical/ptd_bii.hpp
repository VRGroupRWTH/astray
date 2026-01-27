#pragma once

#include <cmath>

#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>

namespace ast::metrics
{
// Implementation based on Motion4D metric by Felix Beslmeisl (tauzero7)
// Reference: https://github.com/tauzero7/Motion4D
//
// Plebanski-Demianski Type B II metric in cylindrical coordinates (t,r,phi,z).
//
// This is an exact solution to Einstein's field equations of Petrov Type D.
// The metric uses hyperbolic coordinates in the r direction with parameter b.
//
// The metric represents a specific class of algebraically special spacetimes.
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class ptd_bii : public metric<coordinate_system_type::cylindrical, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  __device__ termination_reason check_termination(const vector_type& position, [[maybe_unused]] const vector_type& direction) const override
  {
    const auto z = position[3];
    
    // Check for singularity at z = 0
    if (z <= static_cast<scalar_type>(0))
      return termination_reason::spacetime_breakdown;
    
    // Check for singularity at z = b
    if (std::abs(z - b) < static_cast<scalar_type>(1e-10))
      return termination_reason::spacetime_breakdown;
    
    return termination_reason::none;
  }

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto r = position[1];
    const auto z = position[3];
    
    const auto t1 = std::sinh(r);
    const auto t2 = std::cosh(r);
    const auto t4 = b - z;
    const auto t5 = t1 * t1;
    const auto t8 = static_cast<scalar_type>(1) / t1 * t2;
    const auto t9 = static_cast<scalar_type>(1) / z;
    const auto t10 = z * z;
    const auto t19 = static_cast<scalar_type>(1) / t4 * t9 * b / static_cast<scalar_type>(2);

    christoffel_symbols_type symbols;
    symbols.setZero();
    
    // Non-zero Christoffel symbols
    symbols(0, 0, 1) = t1 * t2;
    symbols(0, 0, 3) = t4 * t5;
    symbols(0, 1, 0) = t8;
    symbols(0, 3, 0) = t9;
    
    symbols(1, 0, 0) = t8;
    symbols(1, 1, 3) = -t4;
    symbols(1, 3, 1) = t9;
    
    symbols(2, 2, 3) = t4 / t10 / z * b / static_cast<scalar_type>(2);
    symbols(2, 3, 2) = -t19;
    
    symbols(3, 0, 0) = t9;
    symbols(3, 1, 1) = t9;
    symbols(3, 2, 2) = -t19;
    symbols(3, 3, 3) = t19;
    
    return symbols;
  }

  scalar_type b = static_cast<scalar_type>(1);
};
}
