#pragma once

#include <cmath>

#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>

namespace ast::metrics
{
// Implementation based on Motion4D metric by Thomas Mueller (tauzero7)
// Reference: https://github.com/tauzero7/Motion4D
//
// Straight spinning string metric in cylindrical coordinates (t,ρ,φ,z).
//
// The line element is given by:
//   ds² = -(c dt - a dφ)² + dρ² + k²ρ² dφ² + dz²
//
// with constants a and k > 0.
//
// This represents the spacetime around an infinitely long, straight, spinning cosmic string.
// Reference: Volker Perlick, Living Reviews in Relativity, 7(9), 2004.
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class straight_spinning_string : public metric<coordinate_system_type::cylindrical, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  __device__ termination_reason check_termination(const vector_type& position, [[maybe_unused]] const vector_type& direction) const override
  {
    // Check for singularity at ρ = 0
    if (position[1] <= static_cast<scalar_type>(0))
      return termination_reason::spacetime_breakdown;
    return termination_reason::none;
  }

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto rho = position[1];
    const auto c = consts::speed_of_light;

    const auto t3 = static_cast<scalar_type>(1) / rho;
    const auto t4 = static_cast<scalar_type>(1) / c * a * t3;
    const auto t5 = k * k;

    christoffel_symbols_type symbols;
    symbols.setZero();
    
    // Non-zero Christoffel symbols
    symbols(1, 2, 0) = t4;
    symbols(1, 2, 2) = t3;
    
    symbols(2, 1, 0) = t4;
    symbols(2, 1, 2) = t3;
    symbols(2, 2, 1) = -t5 * rho;
    
    return symbols;
  }

  void set_parameters(scalar_type a_val, scalar_type k_val)
  {
    a = a_val;
    k = k_val;
  }

  scalar_type a = static_cast<scalar_type>(1);
  scalar_type k = static_cast<scalar_type>(0.5);
};
}
