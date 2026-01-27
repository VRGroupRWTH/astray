#pragma once

#include <cmath>

#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>

namespace ast::metrics
{
// Implementation based on Motion4D metric by Felix Beslmeisl (tauzero7)
// Reference: https://github.com/tauzero7/Motion4D
//
// Plebanski-Demianski Type A I metric in spherical coordinates (t,r,theta,phi).
//
// This is an exact solution to Einstein's field equations of Petrov Type D.
// The line element has a simple form with parameter b.
//
// The metric represents a specific class of algebraically special spacetimes.
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class ptd_ai : public metric<coordinate_system_type::spherical, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  __device__ termination_reason check_termination(const vector_type& position, [[maybe_unused]] const vector_type& direction) const override
  {
    const auto r = position[1];
    
    // Check for singularity at r = 0
    if (r <= static_cast<scalar_type>(0))
      return termination_reason::spacetime_breakdown;
    
    // Check for singularity at r = b
    if (std::abs(r - b) < consts::epsilon)
      return termination_reason::spacetime_breakdown;
    
    return termination_reason::none;
  }

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto r = position[1];
    const auto theta = position[2];
    
    const auto t1 = -r + b;
    const auto t2 = r * r;
    const auto t8 = static_cast<scalar_type>(1) / r;
    const auto t12 = t8 / t1 * b / static_cast<scalar_type>(2);
    const auto t13 = std::sin(theta);
    const auto t15 = std::cos(theta);
    const auto t16 = static_cast<scalar_type>(1) / t13 * t15;
    const auto t17 = t13 * t13;

    christoffel_symbols_type symbols;
    symbols.setZero();
    
    // Non-zero Christoffel symbols
    symbols(0, 0, 1) = -t1 / t2 / r * b / static_cast<scalar_type>(2);
    symbols(0, 1, 0) = -t12;
    
    symbols(1, 0, 0) = -t12;
    symbols(1, 1, 1) = t12;
    symbols(1, 2, 2) = t8;
    symbols(1, 3, 3) = t8;
    
    symbols(2, 1, 2) = t8;
    symbols(2, 2, 1) = t1;
    symbols(2, 3, 3) = t16;
    
    symbols(3, 1, 3) = t8;
    symbols(3, 2, 3) = t16;
    symbols(3, 3, 1) = t1 * t17;
    symbols(3, 3, 2) = -t13 * t15;
    
    return symbols;
  }

  scalar_type b = static_cast<scalar_type>(1);
};
}
