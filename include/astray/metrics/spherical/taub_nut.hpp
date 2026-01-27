#pragma once

#include <cmath>

#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>

namespace ast::metrics
{
// Implementation based on Motion4D metric by Thomas Mueller (tauzero7)
// Reference: https://github.com/tauzero7/Motion4D
//
// Taub-NUT metric in Boyer-Lindquist like spherical coordinates (t,r,theta,phi).
//
// The line element is given by:
//   ds² = -Δ/Σ (dt + 2l cos(θ) dφ)² + Σ² (dr²/Δ + dθ² + sin²(θ) dφ²)
//
// where Δ = r² - 2Mr - l² and Σ = r² + l²
//
// The Taub-NUT spacetime represents a rotating solution with a gravitomagnetic
// monopole charge (NUT parameter l).
//
// Reference:
//   Bini et al, "Circular holonomy in the Taub-NUT spacetime",
//   Class. Quantum Grav. 19, 5481 (2002).
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class taub_nut : public metric<coordinate_system_type::spherical, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  __device__ termination_reason check_termination(const vector_type& position, [[maybe_unused]] const vector_type& direction) const override
  {
    const auto r = position[1];
    const auto l_sq = l * l;
    const auto Delta = r * r - static_cast<scalar_type>(2) * mass * r - l_sq;
    
    // Check for horizon (where Delta = 0)
    if (Delta <= static_cast<scalar_type>(0))
      return termination_reason::spacetime_breakdown;
    
    return termination_reason::none;
  }

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto r     = position[1];
    const auto theta = position[2];
    const auto M     = mass;
    
    const auto r_sq = r * r;
    const auto l_sq = l * l;
    const auto Sigma = r_sq + l_sq;
    const auto Delta = r_sq - static_cast<scalar_type>(2) * M * r - l_sq;
    
    const auto st = std::sin(theta);
    const auto ct = std::cos(theta);
    const auto st_sq = st * st;
    const auto ct_sq = ct * ct;
    const auto tan_theta = std::tan(theta);
    
    christoffel_symbols_type symbols;
    symbols.setZero();
    
    // Non-zero Christoffel symbols (directly from Motion4D)
    symbols(0, 0, 1) = (-static_cast<scalar_type>(2) * r * Delta + (static_cast<scalar_type>(-2) * M + static_cast<scalar_type>(2) * r) * Sigma)
        * Delta / (static_cast<scalar_type>(2) * std::pow(Sigma, static_cast<scalar_type>(3)));
    
    symbols(0, 1, 0) = (-static_cast<scalar_type>(2) * r * Delta + (static_cast<scalar_type>(-2) * M + static_cast<scalar_type>(2) * r) * Sigma)
        / (static_cast<scalar_type>(2) * Sigma * Delta);
    
    symbols(0, 2, 0) = -static_cast<scalar_type>(2) * l_sq * Delta / (Sigma * Sigma * tan_theta);
    symbols(0, 2, 3) = l * Delta / (Sigma * Sigma * st);
    
    symbols(0, 3, 1) = l * (-static_cast<scalar_type>(2) * r * Delta + (static_cast<scalar_type>(-2) * M + static_cast<scalar_type>(2) * r) * Sigma)
        * Delta * ct / std::pow(Sigma, static_cast<scalar_type>(3));
    symbols(0, 3, 2) = -l * Delta * st / (Sigma * Sigma);
    
    symbols(1, 0, 0) = symbols(0, 1, 0);
    
    symbols(1, 1, 1) = (static_cast<scalar_type>(2) * r * Delta - (static_cast<scalar_type>(-2) * M + static_cast<scalar_type>(2) * r) * Sigma)
        / (static_cast<scalar_type>(2) * Sigma * Delta);
    
    symbols(1, 2, 2) = r / Sigma;
    
    symbols(1, 3, 0) = l * (-static_cast<scalar_type>(4) * r * Delta + (static_cast<scalar_type>(-2) * M + static_cast<scalar_type>(2) * r) * Sigma) * ct
        / (Sigma * Delta);
    symbols(1, 3, 3) = r / Sigma;
    
    symbols(2, 0, 0) = symbols(0, 2, 0);
    symbols(2, 0, 3) = symbols(0, 2, 3);
    
    symbols(2, 1, 2) = r / Sigma;
    
    symbols(2, 2, 1) = -r * Delta / Sigma;
    
    symbols(2, 3, 0) = l * (-static_cast<scalar_type>(4) * l_sq * Delta * ct_sq + Sigma * Sigma * st_sq - static_cast<scalar_type>(2) * Sigma * Sigma)
        / (Sigma * Sigma * st);
    symbols(2, 3, 3) = (static_cast<scalar_type>(2) * l_sq * Delta + Sigma * Sigma) / (Sigma * Sigma * tan_theta);
    
    symbols(3, 0, 1) = symbols(0, 3, 1);
    symbols(3, 0, 2) = symbols(0, 3, 2);
    
    symbols(3, 1, 0) = symbols(1, 3, 0);
    symbols(3, 1, 3) = symbols(1, 3, 3);
    
    symbols(3, 2, 0) = symbols(2, 3, 0);
    symbols(3, 2, 3) = symbols(2, 3, 3);
    
    symbols(3, 3, 1) = Delta * (-static_cast<scalar_type>(8) * l_sq * r * Delta * ct_sq 
        + static_cast<scalar_type>(4) * l_sq * (static_cast<scalar_type>(-2) * M + static_cast<scalar_type>(2) * r) * Sigma * ct_sq
        + static_cast<scalar_type>(2) * r * Sigma * Sigma * ct_sq - static_cast<scalar_type>(2) * r * Sigma * Sigma)
        / (static_cast<scalar_type>(2) * std::pow(Sigma, static_cast<scalar_type>(3)));
    symbols(3, 3, 2) = (-static_cast<scalar_type>(4) * l_sq * Delta - Sigma * Sigma) * st * ct / (Sigma * Sigma);
    
    return symbols;
  }

  // Mass parameter
  scalar_type mass = static_cast<scalar_type>(1);
  
  // NUT parameter (gravitomagnetic monopole charge)
  scalar_type l = static_cast<scalar_type>(0.5);
};
}
