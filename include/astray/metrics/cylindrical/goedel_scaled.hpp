#pragma once

#include <cmath>

#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>

namespace ast::metrics
{
// Implementation based on Motion4D metric by Frank Grave / Thomas Mueller (tauzero7)
// Reference: https://github.com/tauzero7/Motion4D
//
// Gödel universe in scaled cylindrical coordinates (T, R, φ, Z).
//
// Scaled coordinates where geodesic shape is independent of rG parameter.
// Uses cylindrical coordinate system with the Gödel radius as a scaling parameter.
//
// Parameters:
//   rG - Gödel radius (geodesic scaling parameter)
//   zeta - rotation parameter
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class goedel_scaled : public metric<coordinate_system_type::cylindrical, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto c = consts::speed_of_light;
    const auto R = position[1];

    const auto t1 = R * R;
    const auto t2 = t1 + static_cast<scalar_type>(1);
    const auto t3 = static_cast<scalar_type>(1) / t2;
    const auto t4 = t3 * R;
    const auto t5 = static_cast<scalar_type>(2) * t4;
    const auto t7 = static_cast<scalar_type>(1) / R * t3;
    const auto t8 = std::sqrt(static_cast<scalar_type>(2));
    const auto t10 = t7 * t8 * c;
    const auto t13 = t2 * t8 * c * R;
    const auto t18 = t8 * t1 * R / c * t3;

    christoffel_symbols_type symbols;
    symbols.setZero();
    
    // Non-zero Christoffel symbols
    symbols(0, 1, 0) = t5;
    symbols(0, 1, 2) = -t10;
    
    symbols(0, 2, 1) = t13;
    
    symbols(1, 0, 0) = t5;
    symbols(1, 0, 2) = -t10;
    symbols(1, 1, 1) = -t4;
    symbols(1, 2, 0) = t18;
    symbols(1, 2, 2) = t7;
    
    symbols(2, 0, 1) = t13;
    symbols(2, 1, 0) = t18;
    symbols(2, 1, 2) = t7;
    symbols(2, 2, 1) = t2 * R * (static_cast<scalar_type>(-1) + static_cast<scalar_type>(2) * t1);
    
    return symbols;
  }

  scalar_type rG = static_cast<scalar_type>(1);
  scalar_type zeta = static_cast<scalar_type>(0);
};
}
