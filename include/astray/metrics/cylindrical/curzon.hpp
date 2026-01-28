#pragma once

#include <cmath>

#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>

namespace ast::metrics
{
// Implementation based on Motion4D metric by Thomas Mueller (tauzero7)
// Reference: https://github.com/tauzero7/Motion4D
//
// Curzon metric in cylindrical coordinates (t, r, φ, z).
//
// Static axisymmetric solution describing the exterior field of 
// an infinitely long line mass along the z-axis.
//
// Parameter:
//   mass - mass parameter m
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class curzon : public metric<coordinate_system_type::cylindrical, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto r = position[1];
    const auto z = position[3];
    const auto m = mass;

    const auto t1 = r * r;
    const auto t2 = z * z;
    const auto t3 = t1 + t2;
    const auto t4 = std::sqrt(t3);
    const auto t18 = t4 * t3;
    const auto t19 = static_cast<scalar_type>(1) / t18;
    const auto t20 = m * t19;
    const auto t21 = t20 * r;
    const auto t23 = t20 * z;
    const auto t26 = m * t18;
    const auto t27 = t26 * t1;
    const auto t29 = t1 * t1;
    const auto t30 = t29 * t1;
    const auto t32 = static_cast<scalar_type>(3) * t29 * t2;
    const auto t33 = t2 * t2;
    const auto t35 = static_cast<scalar_type>(3) * t1 * t33;
    const auto t36 = t33 * t2;
    const auto t38 = t3 * t3;
    const auto t40 = static_cast<scalar_type>(1) / t4 / t38;
    const auto t42 = m * r * (-t27 + t26 * t2 + t30 + t32 + t35 + t36) * t40;
    const auto t47 = m * z * (static_cast<scalar_type>(2) * t27 - t30 - t32 - t35 - t36) * t40;
    const auto t50 = t18 - m * t1;
    const auto t52 = static_cast<scalar_type>(1) / r * t50 * t19;
    
    // Additional terms for metric computation
    const auto t7 = std::exp(m / t4);
    const auto t8 = t7 * t7;
    const auto t9 = t8 * t8;
    const auto t11 = m * m;
    const auto t13 = t38;
    const auto t16 = std::exp(t11 * t1 / t13);
    const auto t17 = static_cast<scalar_type>(1) / t9 * t16;

    christoffel_symbols_type symbols;
    symbols.setZero();
    
    // Non-zero Christoffel symbols
    symbols(0, 0, 1) = t17 * t21;
    symbols(0, 0, 3) = t17 * t23;
    
    symbols(0, 1, 0) = t21;
    symbols(0, 3, 0) = t23;
    
    symbols(1, 0, 0) = t21;
    symbols(1, 1, 1) = -t42;
    symbols(1, 1, 3) = -t47;
    symbols(1, 2, 2) = t52;
    symbols(1, 3, 1) = t47;
    symbols(1, 3, 3) = -t42;
    
    symbols(2, 1, 2) = t52;
    symbols(2, 2, 1) = -t16 * r * t50 * t19;
    symbols(2, 2, 3) = t16 * t1 * t23;
    symbols(2, 3, 2) = -t23;
    
    symbols(3, 0, 0) = t23;
    symbols(3, 1, 1) = t47;
    symbols(3, 1, 3) = -t42;
    symbols(3, 2, 2) = -t23;
    symbols(3, 3, 1) = t42;
    symbols(3, 3, 3) = t47;
    
    return symbols;
  }

  scalar_type mass = static_cast<scalar_type>(1);
};
}
