#pragma once

#include <cmath>

#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>

namespace ast::metrics
{
// Implementation based on Motion4D metric by Felix Beslmeisl (tauzero7)
// Reference: https://github.com/tauzero7/Motion4D
//
// Plebanski-Demianski Type C metric in custom coordinates (t,u,x,y).
//
// This is an exact solution to Einstein's field equations of Petrov Type D.
// The metric uses special coordinates with parameters a and b.
//
// Note: In Motion4D this uses custom coordinates, but we treat it as Cartesian
// with coordinates labeled as (t, u, x, y) for implementation purposes.
//
// The metric represents a specific class of algebraically special spacetimes.
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class ptd_c : public metric<coordinate_system_type::cartesian, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  __device__ termination_reason check_termination(const vector_type& position, [[maybe_unused]] const vector_type& direction) const override
  {
    const auto x = position[2];
    const auto y = position[3];
    
    // Check for singularity at x + y = 0
    if (std::abs(x + y) < static_cast<scalar_type>(1e-10))
      return termination_reason::spacetime_breakdown;
    
    return termination_reason::none;
  }

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto x = position[2];
    const auto y = position[3];
    
    const auto t2 = static_cast<scalar_type>(1) / (x + y);
    const auto t3 = x * x;
    const auto t4 = t3 * x;
    const auto t5 = x * a;
    const auto t6 = t4 + t5 + b;
    const auto t7 = t2 * t6;
    const auto t8 = y * y;
    const auto t9 = t8 * y;
    const auto t10 = y * a;
    const auto t11 = t9 + t10 - b;
    const auto t12 = t7 * t11;
    const auto t13 = t2 * t11;
    const auto t15 = static_cast<scalar_type>(3) * t8 * x;
    const auto t16 = static_cast<scalar_type>(2) * b;
    const auto t17 = t15 + t9 + t5 - t10 + t16;
    const auto t20 = static_cast<scalar_type>(1) / t11;
    const auto t21 = t2 * t20;
    const auto t23 = t21 * t17 / static_cast<scalar_type>(2);
    const auto t25 = static_cast<scalar_type>(3) * t3 * y;
    const auto t26 = -t4 + t5 + t16 - t25 - t10;
    const auto t29 = static_cast<scalar_type>(1) / t6;
    const auto t30 = t2 * t29;
    const auto t32 = t30 * t26 / static_cast<scalar_type>(2);

    christoffel_symbols_type symbols;
    symbols.setZero();
    
    // Non-zero Christoffel symbols
    symbols(0, 0, 2) = t12;
    symbols(0, 0, 3) = t13 * t17 / static_cast<scalar_type>(2);
    symbols(0, 2, 0) = -t2;
    symbols(0, 3, 0) = t23;
    
    symbols(1, 1, 2) = t7 * t26 / static_cast<scalar_type>(2);
    symbols(1, 1, 3) = -t12;
    symbols(1, 2, 1) = -t32;
    symbols(1, 3, 1) = -t2;
    
    symbols(2, 0, 0) = -t2;
    symbols(2, 1, 1) = -t32;
    symbols(2, 2, 2) = -t30 * (static_cast<scalar_type>(5) * t4 + static_cast<scalar_type>(3) * t5 + t16 + t25 + t10) / static_cast<scalar_type>(2);
    symbols(2, 2, 3) = -t13 * t29;
    symbols(2, 3, 2) = -t2;
    symbols(2, 3, 3) = -t2;
    
    symbols(3, 0, 0) = t23;
    symbols(3, 1, 1) = -t2;
    symbols(3, 2, 2) = -t2;
    symbols(3, 2, 3) = -t2;
    symbols(3, 3, 2) = -t7 * t20;
    symbols(3, 3, 3) = -t21 * (static_cast<scalar_type>(5) * t9 + static_cast<scalar_type>(3) * t10 - t16 + t15 + t5) / static_cast<scalar_type>(2);
    
    return symbols;
  }

  scalar_type a = static_cast<scalar_type>(1);
  scalar_type b = static_cast<scalar_type>(1);
};
}
