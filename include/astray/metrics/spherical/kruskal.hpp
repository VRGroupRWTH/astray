#pragma once

#include <cmath>
#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>
#include <astray/math/special_functions.hpp>

namespace ast::metrics
{
// Implementation based on Motion4D metric by Thomas Mueller (tauzero7)
// Reference: https://github.com/tauzero7/Motion4D
//
// Kruskal metric in Kruskal-Szekeres coordinates (T,X,theta,phi).
// This is a maximally extended form of the Schwarzschild solution.
// Uses Lambert W function to convert from Kruskal coordinates to Schwarzschild r.
// The coordinate transformation involves: r = 2M(1 + W((X²-T²)/e))

template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class kruskal : public metric<coordinate_system_type::spherical, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  __device__ __host__ inline scalar_type get_r(scalar_type T, scalar_type X) const
  {
    // Calculate r from Kruskal coordinates using Lambert W function
    // r = rs * (1 + W((X²-T²)/e))
    const auto a = (X * X - T * T) / std::exp(static_cast<scalar_type>(1));
    const auto W = math::lambert_w0(a);
    return rs * (W + static_cast<scalar_type>(1));
  }

  __device__ termination_reason check_termination(const vector_type& position, const vector_type& direction) const override
  {
    const auto T = position[0];
    const auto X = position[1];
    
    // Break condition: inside event horizon region
    if (X * X - T * T < static_cast<scalar_type>(-1))
      return termination_reason::spacetime_breakdown;
    
    return termination_reason::none;
  }

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto T = position[0];
    const auto X = position[1];
    const auto r = get_r(T, X);
    const auto theta = position[2];

    const auto exp_neg_r_rs = std::exp(-r / rs);
    const auto rs_over_r2 = rs / (r * r);
    const auto factor1 = T * rs * (r + rs) / (r * r) * exp_neg_r_rs;
    const auto factor2 = -X * rs * (r + rs) / (r * r) * exp_neg_r_rs;
    const auto factor3 = -static_cast<scalar_type>(2) * rs * rs / (r * r) * T * exp_neg_r_rs;
    const auto factor4 = static_cast<scalar_type>(2) * rs * rs / (r * r) * X * exp_neg_r_rs;
    const auto factor5 = -r / (static_cast<scalar_type>(2) * rs);
    
    const auto sin_theta = std::sin(theta);
    const auto cos_theta = std::cos(theta);
    const auto sin2_theta = sin_theta * sin_theta;

    christoffel_symbols_type symbols;
    symbols.setZero();
    
    // Γ^T_TT, Γ^T_TX, Γ^T_θθ, Γ^T_φφ
    symbols(0, 0, 0) = factor1;
    symbols(0, 0, 1) = factor2;
    symbols(0, 1, 0) = factor2;
    symbols(0, 1, 1) = factor1;
    symbols(0, 2, 2) = factor3;
    symbols(0, 3, 3) = factor3;
    
    // Γ^X_TX, Γ^X_XX, Γ^X_θθ, Γ^X_φφ
    symbols(1, 0, 0) = factor2;
    symbols(1, 0, 1) = factor1;
    symbols(1, 1, 0) = factor1;
    symbols(1, 1, 1) = factor2;
    symbols(1, 2, 2) = factor4;
    symbols(1, 3, 3) = factor4;
    
    // Γ^θ_Tθ, Γ^θ_Xθ, Γ^θ_θT, Γ^θ_θX, Γ^θ_φφ
    symbols(2, 0, 2) = factor3;
    symbols(2, 1, 2) = factor4;
    symbols(2, 2, 0) = factor5 * T;
    symbols(2, 2, 1) = factor5 * X;
    symbols(2, 3, 3) = cos_theta / sin_theta;
    
    // Γ^φ_Tφ, Γ^φ_Xφ, Γ^φ_θφ, Γ^φ_φT, Γ^φ_φX, Γ^φ_φθ
    symbols(3, 0, 3) = factor3;
    symbols(3, 1, 3) = factor4;
    symbols(3, 2, 3) = cos_theta / sin_theta;
    symbols(3, 3, 0) = factor5 * T * sin2_theta;
    symbols(3, 3, 1) = factor5 * X * sin2_theta;
    symbols(3, 3, 2) = -sin_theta * cos_theta;
    
    return symbols;
  }
  
  scalar_type mass = static_cast<scalar_type>(1);
  scalar_type rs = static_cast<scalar_type>(2) * mass;
};
}
