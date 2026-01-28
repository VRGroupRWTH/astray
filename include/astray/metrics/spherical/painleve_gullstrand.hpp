#pragma once

#include <cmath>

#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>

namespace ast::metrics
{
// Implementation based on Motion4D metric by Thomas Mueller (tauzero7)
// Reference: https://github.com/tauzero7/Motion4D
//
// Schwarzschild spacetime in Painlevé-Gullstrand coordinates (T,r,theta,phi).
//
// The line element is given by:
//   ds^2 = -c^2 dT^2 + (dr + sqrt(rs/r) c dT)^2 + r^2 (dθ^2 + sin²θ dφ^2)
//
// where rs = 2GM/c^2 is the Schwarzschild radius.
//
// These coordinates are regular at the Schwarzschild radius and represent
// the metric from the viewpoint of a freely falling observer.
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class painleve_gullstrand : public metric<coordinate_system_type::spherical, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  __device__ termination_reason check_termination(const vector_type& position, [[maybe_unused]] const vector_type& direction) const override
  {
    // The Painlevé-Gullstrand coordinates are regular at the Schwarzschild radius
    // but we still need to check for r = 0 singularity
    if (position[1] <= static_cast<scalar_type>(0))
      return termination_reason::spacetime_breakdown;
    return termination_reason::none;
  }

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto r     = position[1];
    const auto theta = position[2];
    const auto rs    = consts::schwarzschild_radius(mass);
    const auto c     = consts::speed_of_light;

    const auto t1  = static_cast<scalar_type>(1) / r;
    const auto t3  = std::sqrt(rs * t1);
    const auto t5  = r * r;
    const auto t7  = rs / t5;
    const auto t9  = c * t3 * t7 / static_cast<scalar_type>(2);
    const auto t10 = r - rs;
    const auto t14 = c * c;
    const auto t18 = t7 / static_cast<scalar_type>(2);
    const auto t19 = static_cast<scalar_type>(1) / c;
    const auto t24 = t19 * t3;
    const auto t26 = std::sin(theta);
    const auto t28 = std::cos(theta);
    const auto t29 = static_cast<scalar_type>(1) / t26 * t28;
    const auto t30 = t26 * t26;

    christoffel_symbols_type symbols;
    symbols.setZero();
    
    // Non-zero Christoffel symbols
    symbols(0, 0, 0) = t9;
    symbols(0, 0, 1) = t10 / t5 / r * t14 * rs / static_cast<scalar_type>(2);
    
    symbols(0, 1, 0) = t18;
    symbols(0, 1, 1) = -t9;
    
    symbols(1, 0, 0) = t18;
    symbols(1, 0, 1) = -t9;
    symbols(1, 1, 0) = t19 / t3 * t7 / static_cast<scalar_type>(2);
    symbols(1, 1, 1) = -t18;
    symbols(1, 2, 2) = t1;
    symbols(1, 3, 3) = t1;
    
    symbols(2, 1, 2) = t1;
    symbols(2, 2, 0) = -t24 * r;
    symbols(2, 2, 1) = -t10;
    symbols(2, 3, 3) = t29;
    
    symbols(3, 1, 3) = t1;
    symbols(3, 2, 3) = t29;
    symbols(3, 3, 0) = -t24 * r * t30;
    symbols(3, 3, 1) = -t10 * t30;
    symbols(3, 3, 2) = -t26 * t28;
    
    return symbols;
  }

  scalar_type mass = static_cast<scalar_type>(1);
};
}
