#pragma once

#include <cmath>
#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>

namespace ast::metrics
{
// Implementation based on Motion4D metric by Thomas Mueller (tauzero7)
// Reference: https://github.com/tauzero7/Motion4D
//
// Hartle-Thorne metric with Gauss-Bonnet corrections in spherical coordinates.
// This represents a slowly rotating neutron star with higher-order corrections.
//
// NOTE: This is a simplified implementation. The full Hartle-Thorne-GB metric requires
// extensive perturbative calculations and auxiliary functions (htt, hrr, hthth, hphph, etc.)
// that depend on the specific stellar model and Gauss-Bonnet coupling.
//
// Reference: Hartle & Thorne, ApJ 153, 807 (1968)
// Gauss-Bonnet corrections: Yagi, Phys. Rev. D 86, 081504 (2012)

template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class hartle_thorne_gb : public metric<coordinate_system_type::spherical, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  __device__ termination_reason check_termination(const vector_type& position, const vector_type& direction) const override
  {
    const auto r = position[1];
    const auto rs = static_cast<scalar_type>(2) * mass;
    
    if (r < rs * static_cast<scalar_type>(1.1))
      return termination_reason::spacetime_breakdown;
    
    return termination_reason::none;
  }

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto r = position[1];
    const auto theta = position[2];
    
    // NOTE: This is a simplified implementation using Schwarzschild as base.
    // The full Hartle-Thorne-GB metric requires:
    // 1. Background Schwarzschild metric
    // 2. Slow rotation corrections (frame-dragging)
    // 3. Oblateness corrections from stellar structure
    // 4. Gauss-Bonnet higher-curvature corrections
    // Each requires solving coupled differential equations
    
    christoffel_symbols_type symbols;
    symbols.setZero();
    
    // Simplified Schwarzschild-like base (placeholder)
    const auto rs = static_cast<scalar_type>(2) * mass;
    const auto inv_r = static_cast<scalar_type>(1) / r;
    const auto f = static_cast<scalar_type>(1) - rs / r;
    
    // Basic Schwarzschild symbols (corrections not implemented)
    symbols(0, 0, 1) = rs / (static_cast<scalar_type>(2) * r * r * f);
    symbols(0, 1, 0) = symbols(0, 0, 1);
    symbols(1, 0, 0) = rs * f / (static_cast<scalar_type>(2) * r * r);
    symbols(1, 1, 1) = -symbols(0, 0, 1);
    symbols(1, 2, 2) = inv_r;
    symbols(1, 3, 3) = inv_r;
    symbols(2, 1, 2) = inv_r;
    symbols(2, 2, 1) = -(r - rs);
    symbols(2, 3, 3) = std::cos(theta) / std::sin(theta);
    symbols(3, 1, 3) = inv_r;
    symbols(3, 2, 3) = std::cos(theta) / std::sin(theta);
    symbols(3, 3, 1) = -(r - rs) * std::sin(theta) * std::sin(theta);
    symbols(3, 3, 2) = -std::sin(theta) * std::cos(theta);
    
    // TODO: Add Hartle-Thorne corrections for rotation and oblateness
    // TODO: Add Gauss-Bonnet higher-curvature corrections
    
    return symbols;
  }
  
  scalar_type mass = static_cast<scalar_type>(1);
  scalar_type angular_momentum = static_cast<scalar_type>(0);  // Slow rotation
  scalar_type eta = static_cast<scalar_type>(0);               // Gauss-Bonnet coupling
};
}
