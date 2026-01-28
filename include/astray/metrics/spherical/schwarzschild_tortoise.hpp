#pragma once

#include <cmath>
#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>
#include <astray/math/special_functions.hpp>

namespace ast::metrics
{
// Implementation based on Schwarzschild tortoise coordinates
// Reference: Standard GR textbooks (Wald, MTW, Carroll)
//
// Tortoise coordinate: r* = r + 2M ln|r/(2M) - 1|
// Inverse relation uses Lambert W function
// This coordinate system is useful for analyzing wave propagation near black holes
// as it makes the lightcone structure more explicit.

template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class schwarzschild_tortoise : public metric<coordinate_system_type::spherical, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  __device__ __host__ inline scalar_type get_r_from_rstar(scalar_type r_star) const
  {
    // Solve r* = r + 2M ln|r/(2M) - 1| for r using Lambert W
    // This is an approximation - the exact inverse requires solving a transcendental equation
    // For large r*: r ≈ r*
    // For r* near 2M: use iterative solution or approximation
    
    // Simple iterative approach (Newton-Raphson)
    scalar_type r = r_star; // Initial guess
    const auto M = mass;
    const auto rs_val = rs;
    
    for (int i = 0; i < 10; ++i)
    {
      const auto f = r + static_cast<scalar_type>(2) * M * std::log(std::abs(r / rs_val - static_cast<scalar_type>(1))) - r_star;
      const auto fprime = static_cast<scalar_type>(1) + static_cast<scalar_type>(2) * M / (r - rs_val);
      r = r - f / fprime;
      if (std::abs(f) < static_cast<scalar_type>(1e-10))
        break;
    }
    
    return r;
  }

  __device__ termination_reason check_termination(const vector_type& position, const vector_type& direction) const override
  {
    const auto r_star = position[1];
    const auto r = get_r_from_rstar(r_star);
    
    if (r < rs * static_cast<scalar_type>(1.01))
      return termination_reason::spacetime_breakdown;
    
    return termination_reason::none;
  }

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto r_star = position[1];
    const auto r = get_r_from_rstar(r_star);
    const auto theta = position[2];
    
    const auto M = mass;
    const auto f = static_cast<scalar_type>(1) - rs / r;
    const auto df_dr = rs / (r * r);
    const auto dr_drstar = f; // dr/dr* = f
    
    const auto sin_theta = std::sin(theta);
    const auto cos_theta = std::cos(theta);
    const auto sin2_theta = sin_theta * sin_theta;
    
    christoffel_symbols_type symbols;
    symbols.setZero();
    
    // Transform Schwarzschild Christoffel symbols to tortoise coordinates
    // The metric in tortoise coordinates is ds² = -f dt² + f⁻¹ dr*² + r² dΩ²
    
    // Γ^t_tr* = (1/2) * (1/f) * df/dr * dr/dr* = (1/2) * df/dr
    const auto factor_t = static_cast<scalar_type>(0.5) * df_dr;
    
    // Γ^r*_tt = (1/2) * f * df/dr * (dr/dr*)
    const auto factor_r = static_cast<scalar_type>(0.5) * f * df_dr * dr_drstar;
    
    // Γ^r*_r*r* = -(1/2) * (df/dr) * (dr/dr*)² 
    const auto factor_rr = -static_cast<scalar_type>(0.5) * df_dr * dr_drstar * dr_drstar / f;
    
    // Γ^r*_θθ = -r * (dr/dr*)
    const auto factor_theta = -r * dr_drstar;
    
    // Γ^θ_r*θ = 1/r * (dr/dr*)⁻¹ = 1/(r * dr_drstar)
    const auto factor_theta_r = static_cast<scalar_type>(1) / (r * dr_drstar);
    
    symbols(0, 0, 1) = factor_t;
    symbols(0, 1, 0) = factor_t;
    
    symbols(1, 0, 0) = factor_r;
    symbols(1, 1, 1) = factor_rr;
    symbols(1, 2, 2) = factor_theta;
    symbols(1, 3, 3) = factor_theta * sin2_theta;
    
    symbols(2, 1, 2) = factor_theta_r;
    symbols(2, 2, 1) = factor_theta_r;
    symbols(2, 3, 3) = cos_theta / sin_theta;
    
    symbols(3, 1, 3) = factor_theta_r;
    symbols(3, 3, 1) = factor_theta_r;
    symbols(3, 3, 2) = -sin_theta * cos_theta;
    
    return symbols;
  }
  
  scalar_type mass = static_cast<scalar_type>(1);
  scalar_type rs = static_cast<scalar_type>(2) * mass;
};
}
