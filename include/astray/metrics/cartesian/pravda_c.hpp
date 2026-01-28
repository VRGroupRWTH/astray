#pragma once

#include <cmath>

#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>

namespace ast::metrics
{
// Implementation based on Motion4D metric by Felix Beslmeisl and Thomas Mueller
// Reference: https://github.com/tauzero7/Motion4D
//
// Pravda C-metric in Cartesian-like coordinates (t,x,y,phi).
//
// This metric represents a pair of accelerating black holes.
// The line element is given in terms of coordinates x and y with 
// parameters A (acceleration) and m (mass).
//
// The C-metric is an important exact solution in general relativity
// describing uniformly accelerating black holes.
//
// Note: Uses custom coordinates (t, x, y, phi) rather than standard Cartesian.
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class pravda_c : public metric<coordinate_system_type::cartesian, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto x = position[1];
    const auto y = position[2];
    const auto A = acceleration;
    const auto m = mass;
    
    const auto mA = m * A;
    
    const auto t2 = static_cast<scalar_type>(1) / (x + y);
    const auto t3 = x * x;
    const auto t4 = mA;
    const auto t6 = t4 * t3 * x;
    const auto t8 = static_cast<scalar_type>(-1) + t3 + static_cast<scalar_type>(2) * t6;
    const auto t9 = t2 * t8;
    const auto t10 = y * y;
    const auto t12 = t4 * t10 * y;
    const auto t14 = static_cast<scalar_type>(1) - t10 + static_cast<scalar_type>(2) * t12;
    const auto t15 = t9 * t14;
    const auto t16 = t2 * t14;
    const auto t17 = y * x;
    const auto t20 = static_cast<scalar_type>(3) * t4 * t10 * x;
    const auto t21 = static_cast<scalar_type>(-1) + t12 - t17 + t20;
    const auto t23 = static_cast<scalar_type>(1) / t14;
    const auto t24 = t2 * t23;
    const auto t25 = t24 * t21;
    const auto t26 = static_cast<scalar_type>(1) / t8;
    const auto t27 = t2 * t26;
    const auto t32 = static_cast<scalar_type>(3) * t4 * y * t3;
    const auto t36 = static_cast<scalar_type>(1) + t6 + t17 + t32;
    const auto t37 = t27 * t36;

    christoffel_symbols_type symbols;
    symbols.setZero();
    
    // Non-zero Christoffel symbols
    symbols(0, 0, 1) = -t15;
    symbols(0, 0, 2) = t16 * t21;
    
    symbols(0, 1, 0) = -t2;
    
    symbols(0, 2, 0) = t25;
    
    symbols(1, 0, 0) = -t2;
    
    symbols(1, 1, 1) = -t27 * (static_cast<scalar_type>(-1) + static_cast<scalar_type>(2) * t3 + static_cast<scalar_type>(5) * t6 + t17 + t32);
    symbols(1, 1, 2) = t16 * t26;
    
    symbols(1, 2, 1) = -t2;
    symbols(1, 2, 2) = -t2;
    
    symbols(1, 3, 3) = t37;
    
    symbols(2, 0, 0) = t25;
    
    symbols(2, 1, 1) = -t2;
    symbols(2, 1, 2) = -t2;
    
    symbols(2, 2, 1) = t9 * t23;
    symbols(2, 2, 2) = -t24 * (static_cast<scalar_type>(1) - static_cast<scalar_type>(2) * t10 + static_cast<scalar_type>(5) * t12 - t17 + t20);
    
    symbols(2, 3, 3) = -t2;
    
    symbols(3, 1, 3) = t37;
    
    symbols(3, 2, 3) = -t2;
    
    symbols(3, 3, 1) = -t9 * t36;
    symbols(3, 3, 2) = t15;
    
    return symbols;
  }

  // Acceleration parameter
  scalar_type acceleration = static_cast<scalar_type>(1);
  
  // Mass parameter
  scalar_type mass = static_cast<scalar_type>(1);
};
}
