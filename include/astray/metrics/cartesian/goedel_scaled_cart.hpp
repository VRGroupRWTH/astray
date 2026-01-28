#pragma once

#include <cmath>

#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>

namespace ast::metrics
{
// Implementation based on Motion4D metric by Frank Grave / Thomas Mueller (tauzero7)
// Reference: https://github.com/tauzero7/Motion4D
//
// Gödel universe in scaled Cartesian coordinates (T, X, Y, Z).
//
// Scaled Cartesian form where geodesic shape is independent of rG parameter.
// Valid at X=Y=0 (no coordinate singularity).
//
// Parameters:
//   rG - Gödel radius (geodesic scaling parameter)
//   zeta - rotation parameter
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class goedel_scaled_cart : public metric<coordinate_system_type::cartesian, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto c = consts::speed_of_light;
    const auto X = position[1];
    const auto Y = position[2];

    const auto t1 = X * X;
    const auto t2 = Y * Y;
    const auto t4 = static_cast<scalar_type>(1) / (static_cast<scalar_type>(1) + t1 + t2);
    const auto t5 = X * t4;
    const auto t6 = static_cast<scalar_type>(2) * t5;
    const auto t10 = std::sqrt(static_cast<scalar_type>(2));
    const auto t13 = (static_cast<scalar_type>(2) + t1 + t2) * X * Y * t4 * t10 * c;
    const auto t14 = t2 * t1;
    const auto t15 = t2 * t2;
    const auto t19 = t10 * c;
    const auto t20 = (t14 + static_cast<scalar_type>(1) + t15 + static_cast<scalar_type>(2) * t2) * t4 * t19;
    const auto t22 = static_cast<scalar_type>(2) * Y * t4;
    const auto t23 = t1 * t1;
    const auto t27 = (t23 + t14 + static_cast<scalar_type>(1) + static_cast<scalar_type>(2) * t1) * t4 * t19;
    const auto t29 = static_cast<scalar_type>(1) / c;
    const auto t32 = static_cast<scalar_type>(2) * Y * t10 * t5 * t29;
    const auto t33 = static_cast<scalar_type>(2) * t14;
    const auto t34 = static_cast<scalar_type>(3) * t2;
    const auto t35 = static_cast<scalar_type>(2) * t15;
    const auto t39 = t33 + t35 + static_cast<scalar_type>(1) + t34;
    const auto t45 = (t1 - t2) * t10 * t4 * t29;
    const auto t46 = static_cast<scalar_type>(2) * t23;
    const auto t47 = static_cast<scalar_type>(3) * t1;
    const auto t48 = t46 + t47 + t33 + static_cast<scalar_type>(1);
    const auto t50 = Y * t48 * t4;
    const auto t52 = X * t39 * t4;

    christoffel_symbols_type symbols;
    symbols.setZero();
    
    // Non-zero Christoffel symbols
    symbols(0, 1, 0) = t6;
    symbols(0, 1, 1) = -t13;
    symbols(0, 1, 2) = -t20;
    
    symbols(0, 2, 0) = t22;
    symbols(0, 2, 1) = t27;
    symbols(0, 2, 2) = t13;
    
    symbols(1, 0, 0) = t6;
    symbols(1, 0, 1) = -t13;
    symbols(1, 0, 2) = -t20;
    symbols(1, 1, 0) = -t32;
    symbols(1, 1, 1) = X * (t33 + t34 + t35 - static_cast<scalar_type>(1)) * t4;
    symbols(1, 1, 2) = Y * t39 * t4;
    symbols(1, 2, 0) = t45;
    symbols(1, 2, 1) = -t50;
    symbols(1, 2, 2) = -t52;
    
    symbols(2, 0, 0) = t22;
    symbols(2, 0, 1) = t27;
    symbols(2, 0, 2) = t13;
    symbols(2, 1, 0) = t45;
    symbols(2, 1, 1) = -t50;
    symbols(2, 1, 2) = -t52;
    symbols(2, 2, 0) = t32;
    symbols(2, 2, 1) = X * t48 * t4;
    symbols(2, 2, 2) = Y * (t46 + t47 + t33 - static_cast<scalar_type>(1)) * t4;
    
    return symbols;
  }

  scalar_type rG = static_cast<scalar_type>(1);
  scalar_type zeta = static_cast<scalar_type>(0);
};
}
