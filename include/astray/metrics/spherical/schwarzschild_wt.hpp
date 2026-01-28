#pragma once

#include <cmath>
#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>

namespace ast::metrics
{
// Implementation based on Motion4D metric by Thomas Mueller (tauzero7)
// Reference: https://github.com/tauzero7/Motion4D
//
// Schwarzschild-WT (Wheeler-Thorne) metric in spherical coordinates (t,r,theta,phi).
//
// Line element:
//   ds^2 = -c^2 dt^2 + dr^2/(1-rs/r) + r^2(dθ^2 + sin²θ dφ^2)
//
// where rs = 2GM/c^2 is the Schwarzschild radius.
//
// This differs from the standard Schwarzschild metric by having a constant
// g_tt component (no (1-rs/r) factor), which makes it useful for Wheeler-Thorne
// coordinates in the study of black hole physics.

template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class schwarzschild_wt : public metric<coordinate_system_type::spherical, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  __device__ termination_reason check_termination(const vector_type& position, const vector_type& direction) const override
  {
    const auto rs = consts::schwarzschild_radius(mass);
    if (position[1] < static_cast<scalar_type>(0) || 
        static_cast<scalar_type>(std::pow(position[1], 2)) <= (static_cast<scalar_type>(1) + consts::epsilon) * static_cast<scalar_type>(std::pow(rs, 2)))
      return termination_reason::spacetime_breakdown;
    return termination_reason::none;
  }

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto rs    = consts::schwarzschild_radius(mass);
    const auto r     = position[1];
    const auto theta = position[2];
    const auto t1    = r - rs;
    const auto t10   = static_cast<scalar_type>(1) / r;
    const auto t14   = t10 / t1 * rs * static_cast<scalar_type>(0.5);
    const auto t15   = std::sin(theta);
    const auto t17   = std::cos(theta);
    const auto t18   = static_cast<scalar_type>(1) / t15 * t17;
    const auto t19   = static_cast<scalar_type>(std::pow(t15, 2));

    christoffel_symbols_type symbols;
    symbols.setZero();
    
    // Note: In Wheeler-Thorne coordinates, g_tt = -c^2 (constant), 
    // so christoffel[0][0][1], christoffel[0][1][0], christoffel[1][0][0] = 0
    // (unlike standard Schwarzschild where they are non-zero)
    symbols(1, 1, 1) = -t14;
    symbols(1, 2, 2) =  t10;
    symbols(1, 3, 3) =  t10;
    symbols(2, 1, 2) =  t10;
    symbols(2, 2, 1) = -t1;
    symbols(2, 3, 3) =  t18;
    symbols(3, 1, 3) =  t10;
    symbols(3, 2, 3) =  t18;
    symbols(3, 3, 1) = -t1  * t19;
    symbols(3, 3, 2) = -t15 * t17;
    return symbols;
  }
  
  scalar_type mass = static_cast<scalar_type>(1);
};
}
