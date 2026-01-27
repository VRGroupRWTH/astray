#pragma once

#include <cmath>

#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>

namespace ast::metrics
{
// Implementation based on Motion4D metric by Frank Grave / Thomas Mueller (tauzero7)
// Reference: https://github.com/tauzero7/Motion4D
//
// Gödel universe in Cartesian coordinates (t, x, y, z).
//
// The metric describes a rotating cosmological solution that allows for closed
// timelike curves. It has cylindrical symmetry around the z-axis.
//
// Parameters:
//   a - Gödel radius parameter (rG = 2a)
//   zeta - rotation parameter
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class goedel_cart : public metric<coordinate_system_type::cartesian, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto c = consts::speed_of_light;
    const auto x = position[1];
    const auto y = position[2];

    const auto t1 = a * a;
    const auto t3 = x * x;
    const auto t4 = y * y;
    const auto t6 = static_cast<scalar_type>(1) / (static_cast<scalar_type>(4) * t1 + t3 + t4);
    const auto t8 = static_cast<scalar_type>(2) * x * t6;
    const auto t14 = static_cast<scalar_type>(1) / t1 / a;
    const auto t16 = std::sqrt(static_cast<scalar_type>(2));
    const auto t20 = y * x * (static_cast<scalar_type>(8) * t1 + t3 + t4) * t14 * t6 * c * t16 / static_cast<scalar_type>(8);
    const auto t21 = t3 * t4;
    const auto t22 = t4 * t4;
    const auto t23 = t1 * t1;
    const auto t24 = static_cast<scalar_type>(16) * t23;
    const auto t25 = t1 * t4;
    const auto t30 = t6 * c * t16;
    const auto t32 = (t21 + t22 + t24 + static_cast<scalar_type>(8) * t25) * t14 * t30 / static_cast<scalar_type>(8);
    const auto t34 = static_cast<scalar_type>(2) * y * t6;
    const auto t35 = t3 * t3;
    const auto t36 = t1 * t3;
    const auto t41 = (t35 + t21 + t24 + static_cast<scalar_type>(8) * t36) * t14 * t30 / static_cast<scalar_type>(8);
    const auto t47 = t6 / a / c;
    const auto t48 = t16 * y * x * t47;
    const auto t49 = static_cast<scalar_type>(8) * t23;
    const auto t50 = static_cast<scalar_type>(6) * t25;
    const auto t54 = static_cast<scalar_type>(1) / t23 * t6;
    const auto t57 = t49 + t50 + t21 + t22;
    const auto t64 = (t3 - t4) * t16 * t47 / static_cast<scalar_type>(2);
    const auto t65 = static_cast<scalar_type>(6) * t36;
    const auto t66 = t49 + t65 + t35 + t21;
    const auto t69 = t66 * y * t54 / static_cast<scalar_type>(8);
    const auto t72 = t57 * x * t54 / static_cast<scalar_type>(8);

    christoffel_symbols_type symbols;
    symbols.setZero();
    
    // Non-zero Christoffel symbols
    symbols(0, 1, 0) = t8;
    symbols(0, 1, 1) = -t20;
    symbols(0, 1, 2) = -t32;
    
    symbols(0, 2, 0) = t34;
    symbols(0, 2, 1) = t41;
    symbols(0, 2, 2) = t20;
    
    symbols(1, 0, 0) = t8;
    symbols(1, 0, 1) = -t20;
    symbols(1, 0, 2) = -t32;
    symbols(1, 1, 0) = -t48;
    symbols(1, 1, 1) = -(t49 - t50 - t21 - t22) * x * t54 / static_cast<scalar_type>(8);
    symbols(1, 1, 2) = t57 * y * t54 / static_cast<scalar_type>(8);
    symbols(1, 2, 0) = t64;
    symbols(1, 2, 1) = -t69;
    symbols(1, 2, 2) = -t72;
    
    symbols(2, 0, 0) = t34;
    symbols(2, 0, 1) = t41;
    symbols(2, 0, 2) = t20;
    symbols(2, 1, 0) = t64;
    symbols(2, 1, 1) = -t69;
    symbols(2, 1, 2) = -t72;
    symbols(2, 2, 0) = t48;
    symbols(2, 2, 1) = t66 * x * t54 / static_cast<scalar_type>(8);
    symbols(2, 2, 2) = -(t49 - t65 - t35 - t21) * y * t54 / static_cast<scalar_type>(8);
    
    return symbols;
  }

  scalar_type a = static_cast<scalar_type>(1);
  scalar_type zeta = static_cast<scalar_type>(0);
};
}
