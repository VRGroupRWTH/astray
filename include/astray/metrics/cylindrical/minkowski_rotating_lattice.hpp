#pragma once

#include <cmath>

#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>

namespace ast::metrics
{
// Implementation based on Motion4D metric by Thomas Mueller (tauzero7)
// Reference: https://github.com/tauzero7/Motion4D
//
// Minkowski spacetime in rotating lattice cylindrical coordinates (t, r, φ, z).
//
// This represents flat Minkowski spacetime from the perspective of a rotating 
// cylindrical coordinate system with angular velocity ω.
//
// Parameter:
//   omega - angular velocity of rotation
//
// Note: Only valid for r < c/ω (outside this region, coordinates become spacelike)
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class minkowski_rotating_lattice : public metric<coordinate_system_type::cylindrical, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  __device__ termination_reason check_termination(const vector_type& position, [[maybe_unused]] const vector_type& direction) const override
  {
    const auto c = consts::speed_of_light;
    const auto r = position[1];
    
    // Coordinates become invalid at r >= c/ω
    if (r >= std::abs(c / omega) - consts::epsilon)
      return termination_reason::spacetime_breakdown;
    return termination_reason::none;
  }

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto r = position[1];

    const auto t1 = omega * omega;
    const auto t3 = static_cast<scalar_type>(1) / r;
    const auto t4 = omega * t3;
    const auto t5 = omega * r;

    christoffel_symbols_type symbols;
    symbols.setZero();
    
    // Non-zero Christoffel symbols
    symbols(0, 0, 1) = -r * t1;
    
    symbols(0, 1, 2) = t4;
    
    symbols(0, 2, 1) = -t5;
    
    symbols(1, 0, 2) = t4;
    symbols(1, 2, 2) = t3;
    
    symbols(2, 0, 1) = -t5;
    symbols(2, 1, 2) = t3;
    symbols(2, 2, 1) = -r;
    
    return symbols;
  }

  scalar_type omega = static_cast<scalar_type>(1);
};
}
