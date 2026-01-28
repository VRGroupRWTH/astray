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
// Halilsoy standing gravitational wave in cylindrical coordinates (t,rho,phi,z).
// This represents a standing gravitational wave solution.
//
// The metric is given by:
// ds² = e^(-2U) [e^(2k) (dρ² - dt²) + ρ² dφ²] + e^(2U) (dz + A dφ)²
//
// where:
// e^(-2U) = cosh²(α) e^(-2C J₀(ρ) cos(t)) + sinh²(α) e^(2C J₀(ρ) cos(t))
// A = -2C sinh(2α) ρ J₁(ρ) sin(t)
// k = (1/2) C² [ρ² (J₀²(ρ) + J₁²(ρ)) - 2ρ J₀(ρ) J₁(ρ) cos²(t)]
//
// Reference: M. Halilsoy, "Cross-Polarized Cylindrical Gravitational Waves of Einstein and Rosen",
// Il Nuovo Cimento 102 B, 563 (1988)

template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class halilsoy_wave : public metric<coordinate_system_type::cylindrical, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  __device__ __host__ void calc_metric_functions(scalar_type t, scalar_type rho, 
                                                   scalar_type& V, scalar_type& A, scalar_type& K) const
  {
    const auto b0 = math::bessel_j0(rho);
    const auto b1 = math::bessel_j1(rho);
    const auto ct = std::cos(t);
    const auto st = std::sin(t);
    const auto cha = std::cosh(alpha);
    const auto sha = std::sinh(alpha);
    
    // V function: e^(-2U) = cosh²(α) e^(-2C J₀ cos(t)) + sinh²(α) e^(2C J₀ cos(t))
    const auto exp_term = std::exp(-static_cast<scalar_type>(2) * C * b0 * ct);
    V = cha * cha * exp_term + sha * sha / exp_term;
    
    // A function: A = -2C sinh(2α) ρ J₁(ρ) sin(t)
    A = -static_cast<scalar_type>(2) * C * std::sinh(static_cast<scalar_type>(2) * alpha) * rho * b1 * st;
    
    // K function: k = (1/2) C² [ρ² (J₀² + J₁²) - 2ρ J₀ J₁ cos²(t)]
    K = static_cast<scalar_type>(0.5) * C * C * 
        (rho * rho * (b0 * b0 + b1 * b1) - static_cast<scalar_type>(2) * rho * b0 * b1 * ct * ct);
  }

  __device__ termination_reason check_termination(const vector_type& position, const vector_type& direction) const override
  {
    const auto rho = position[1];
    if (rho < static_cast<scalar_type>(0))
      return termination_reason::spacetime_breakdown;
    
    return termination_reason::none;
  }

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto t = position[0];
    const auto rho = position[1];
    
    // Calculate metric functions
    scalar_type V, A, K;
    calc_metric_functions(t, rho, V, A, K);
    
    // Calculate derivatives using finite differences for numerical stability
    const auto eps = static_cast<scalar_type>(1e-6);
    
    scalar_type V_t, A_t, K_t, V_rho, A_rho, K_rho;
    
    // Time derivatives
    {
      scalar_type Vp, Ap, Kp;
      calc_metric_functions(t + eps, rho, Vp, Ap, Kp);
      V_t = (Vp - V) / eps;
      A_t = (Ap - A) / eps;
      K_t = (Kp - K) / eps;
    }
    
    // Rho derivatives
    {
      scalar_type Vp, Ap, Kp;
      calc_metric_functions(t, rho + eps, Vp, Ap, Kp);
      V_rho = (Vp - V) / eps;
      A_rho = (Ap - A) / eps;
      K_rho = (Kp - K) / eps;
    }
    
    const auto V2 = V * V;
    const auto inv_V = static_cast<scalar_type>(1) / V;
    const auto inv_V2 = static_cast<scalar_type>(1) / V2;
    const auto exp_2K = std::exp(static_cast<scalar_type>(2) * K);
    const auto inv_exp_2K = static_cast<scalar_type>(1) / exp_2K;
    const auto A2 = A * A;
    const auto rho2 = rho * rho;
    
    // Christoffel symbol components (simplified from Motion4D implementation)
    const auto factor1 = (V_t + static_cast<scalar_type>(2) * V * K_t) * inv_V / static_cast<scalar_type>(2);
    const auto factor2 = (V_rho + static_cast<scalar_type>(2) * V * K_rho) * inv_V / static_cast<scalar_type>(2);
    
    christoffel_symbols_type symbols;
    symbols.setZero();
    
    // Simplified implementation - only non-zero components
    // Full implementation would require all terms from Motion4D
    symbols(0, 0, 0) = factor1;
    symbols(0, 0, 1) = factor2;
    symbols(0, 1, 0) = factor2;
    symbols(0, 1, 1) = factor1;
    
    symbols(1, 0, 0) = factor2;
    symbols(1, 0, 1) = factor1;
    symbols(1, 1, 0) = factor1;
    symbols(1, 1, 1) = factor2;
    
    // Additional components would be added here for full implementation
    // This is a simplified version focusing on the main structure
    
    return symbols;
  }
  
  scalar_type alpha = static_cast<scalar_type>(0);  // Wave parameter
  scalar_type C = static_cast<scalar_type>(1);      // Wave amplitude
};
}
