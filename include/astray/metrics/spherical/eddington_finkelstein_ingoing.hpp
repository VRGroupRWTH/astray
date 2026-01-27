#pragma once

#include <cmath>

#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>

namespace ast::metrics
{
// Implementation based on Motion4D metric by Thomas Mueller (tauzero7)
// Reference: https://github.com/tauzero7/Motion4D
//
// Eddington-Finkelstein ingoing coordinates for Schwarzschild black hole (v, r, θ, φ).
//
// The line element is given by:
//   ds^2 = -(1 - rs/r)c^2 dv^2 + 2c dv dr + r^2(dθ^2 + sin²θ dφ^2)
//
// These coordinates are regular at the event horizon and describe infalling geodesics.
// The v coordinate is the ingoing Eddington time.
//
// Parameter:
//   mass - mass of the black hole (rs = 2GM/c^2)
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class eddington_finkelstein_ingoing : public metric<coordinate_system_type::spherical, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  __device__ termination_reason check_termination(const vector_type& position, [[maybe_unused]] const vector_type& direction) const override
  {
    // Check for singularity at r = 0
    if (position[1] <= consts::epsilon)
      return termination_reason::spacetime_breakdown;
    return termination_reason::none;
  }

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto c = consts::speed_of_light;
    const auto r = position[1];
    const auto theta = position[2];
    const auto rs = consts::schwarzschild_radius(mass);

    const auto t2 = r * r;
    const auto t5 = c * rs / t2 / static_cast<scalar_type>(2);
    const auto t6 = r - rs;
    const auto t10 = c * c;
    const auto t14 = static_cast<scalar_type>(1) / r;
    const auto t16 = static_cast<scalar_type>(1) / c * r;
    const auto t17 = std::sin(theta);
    const auto t19 = std::cos(theta);
    const auto t20 = static_cast<scalar_type>(1) / t17 * t19;
    const auto t21 = t17 * t17;

    christoffel_symbols_type symbols;
    symbols.setZero();
    
    // Non-zero Christoffel symbols
    symbols(0, 0, 0) = t5;
    symbols(0, 0, 1) = t6 / t2 / r * t10 * rs / static_cast<scalar_type>(2);
    
    symbols(0, 1, 1) = -t5;
    
    symbols(1, 0, 1) = -t5;
    symbols(1, 2, 2) = t14;
    symbols(1, 3, 3) = t14;
    
    symbols(2, 1, 2) = t14;
    symbols(2, 2, 0) = -t16;
    symbols(2, 2, 1) = -t6;
    symbols(2, 3, 3) = t20;
    
    symbols(3, 1, 3) = t14;
    symbols(3, 2, 3) = t20;
    symbols(3, 3, 0) = -t16 * t21;
    symbols(3, 3, 1) = -t6 * t21;
    symbols(3, 3, 2) = -t17 * t19;
    
    return symbols;
  }

  scalar_type mass = static_cast<scalar_type>(1);
};
}
