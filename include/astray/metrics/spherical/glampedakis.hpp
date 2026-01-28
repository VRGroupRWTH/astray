#pragma once

#include <cmath>
#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>

namespace ast::metrics
{
// Implementation based on Motion4D metric by Thomas Mueller (tauzero7)
// Reference: https://github.com/tauzero7/Motion4D
//
// Glampedakis metric in spherical coordinates.
// This represents a slowly rotating neutron star with perturbative corrections
// beyond the Kerr metric, accounting for relativistic effects in the stellar structure.
//
// NOTE: This is a simplified implementation. The full Glampedakis metric requires:
// 1. Background Kerr metric components
// 2. Perturbative corrections from stellar structure (similar to Hartle-Thorne)
// 3. Complex auxiliary functions and their derivatives
// 4. Numerical integration of stellar structure equations
//
// A complete implementation would require extensive numerical calculations that
// are beyond the scope of a basic implementation.
//
// Reference: Glampedakis & Babak, Classical and Quantum Gravity 23, 4167 (2006)

template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class glampedakis : public metric<coordinate_system_type::spherical, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  __device__ termination_reason check_termination(const vector_type& position, const vector_type& direction) const override
  {
    const auto r = position[1];
    const auto rs = static_cast<scalar_type>(2) * mass;
    
    // Simple horizon check (would need refinement for rotating case)
    if (r < rs * static_cast<scalar_type>(1.1))
      return termination_reason::spacetime_breakdown;
    
    return termination_reason::none;
  }

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto r = position[1];
    const auto theta = position[2];
    
    // NOTE: This is a simplified Kerr-like implementation as a placeholder.
    // The full Glampedakis metric requires:
    // 1. Kerr background metric g^Kerr_μν(r, θ, a)
    // 2. Perturbative corrections h_μν(r, θ, ε) from stellar structure
    // 3. Combined metric g_μν = g^Kerr_μν + ε h_μν
    // 4. Christoffel symbols computed from the combined metric
    //
    // Each perturbative correction involves solving coupled PDEs for the
    // interior and exterior of the star.
    
    christoffel_symbols_type symbols;
    symbols.setZero();
    
    // Simplified Kerr-like base (Schwarzschild limit for a=0)
    const auto rs = static_cast<scalar_type>(2) * mass;
    const auto a = angular_momentum;
    const auto a2 = a * a;
    const auto r2 = r * r;
    const auto cos_theta = std::cos(theta);
    const auto sin_theta = std::sin(theta);
    const auto sin2_theta = sin_theta * sin_theta;
    
    // Simplified expressions (Schwarzschild-like)
    const auto inv_r = static_cast<scalar_type>(1) / r;
    const auto f = static_cast<scalar_type>(1) - rs / r;
    
    // Basic Schwarzschild symbols as placeholder
    symbols(0, 0, 1) = rs / (static_cast<scalar_type>(2) * r2 * f);
    symbols(0, 1, 0) = symbols(0, 0, 1);
    
    symbols(1, 0, 0) = rs * f / (static_cast<scalar_type>(2) * r2);
    symbols(1, 1, 1) = -symbols(0, 0, 1);
    symbols(1, 2, 2) = inv_r;
    symbols(1, 3, 3) = inv_r;
    
    symbols(2, 1, 2) = inv_r;
    symbols(2, 2, 1) = -(r - rs);
    symbols(2, 3, 3) = cos_theta / sin_theta;
    
    symbols(3, 1, 3) = inv_r;
    symbols(3, 2, 3) = cos_theta / sin_theta;
    symbols(3, 3, 1) = -(r - rs) * sin2_theta;
    symbols(3, 3, 2) = -sin_theta * cos_theta;
    
    // TODO: Add Kerr corrections for rotation (g_tφ terms)
    // TODO: Add Glampedakis perturbative corrections for stellar structure
    // These would involve:
    // - Frame dragging effects from rotation
    // - Multipole moments of the stellar matter distribution
    // - Relativistic corrections to the stellar structure
    
    return symbols;
  }
  
  scalar_type mass = static_cast<scalar_type>(1);
  scalar_type angular_momentum = static_cast<scalar_type>(0);  // Slow rotation parameter
  scalar_type epsilon = static_cast<scalar_type>(0);           // Perturbation parameter
};
}
