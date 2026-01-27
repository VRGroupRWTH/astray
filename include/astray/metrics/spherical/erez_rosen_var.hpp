#pragma once

#include <cmath>

#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>

namespace ast::metrics
{
// Implementation based on Motion4D metric by Thomas Mueller (tauzero7)
// Reference: https://github.com/tauzero7/Motion4D
//
// Erez-Rosen metric with variable deformation parameter in spherical coordinates (t,r,theta,phi).
//
// This is a static axisymmetric solution representing a deformed mass distribution.
// The metric includes potentials ψ, γ, and Δ that depend on position:
//
//   Δ = r² - 2Mr + M²sin²θ
//   γ = (1/2)log[(r² - 2Mr)/Δ]
//   ψ = (1/2)log(1 - 2M/r)
//
// The line element involves these potentials combined with standard spherical terms.
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class erez_rosen_var : public metric<coordinate_system_type::spherical, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  __device__ termination_reason check_termination(const vector_type& position, [[maybe_unused]] const vector_type& direction) const override
  {
    const auto r = position[1];
    
    // Check for singularities
    if (r <= static_cast<scalar_type>(2) * mass)
      return termination_reason::spacetime_breakdown;
    
    return termination_reason::none;
  }

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto r     = position[1];
    const auto theta = position[2];
    const auto m     = mass;
    
    // Calculate potentials
    const auto st = std::sin(theta);
    const auto ct = std::cos(theta);
    const auto st_sq = st * st;
    
    const auto Delta = r * r - static_cast<scalar_type>(2) * m * r + m * m * st_sq;
    const auto g = static_cast<scalar_type>(0.5) * std::log((r * r - static_cast<scalar_type>(2) * m * r) / Delta);
    const auto psi = static_cast<scalar_type>(0.5) * std::log(static_cast<scalar_type>(1) - static_cast<scalar_type>(2) * m / r);
    
    // Calculate derivatives
    const auto dDdr = static_cast<scalar_type>(2) * r - static_cast<scalar_type>(2) * m;
    const auto dDdtheta = static_cast<scalar_type>(2) * m * m * st * ct;
    
    const auto dgdr = ((m * m * r - m * m * m) * st_sq)
        / ((m * m * r * r - static_cast<scalar_type>(2) * m * m * m * r) * st_sq + r * r * r * r 
           - static_cast<scalar_type>(4) * m * r * r * r + static_cast<scalar_type>(4) * m * m * r * r);
    const auto dgdtheta = -(m * m * ct * st) / (m * m * st_sq + r * r - static_cast<scalar_type>(2) * m * r);
    
    const auto dpsidr = m / (r * r - static_cast<scalar_type>(2) * m * r);
    const auto dpsidtheta = static_cast<scalar_type>(0);
    
    // Precompute common terms
    const auto t2 = -r + static_cast<scalar_type>(2) * m;
    const auto t3 = r * t2;
    const auto t4 = psi;
    const auto t5 = std::exp(t4);
    const auto t6 = t5 * t5;
    const auto t7 = t6 * t6;
    const auto t9 = g;
    const auto t10 = std::exp(t9);
    const auto t11 = t10 * t10;
    const auto t12 = static_cast<scalar_type>(1) / t11;
    const auto t13 = Delta;
    const auto t14 = static_cast<scalar_type>(1) / t13;
    const auto t15 = t12 * t14;
    const auto t16 = dpsidr;
    const auto t20 = dpsidtheta;
    const auto t25 = static_cast<scalar_type>(1) / r / t2;
    const auto t26 = dgdr;
    const auto t27 = t13 * t26;
    const auto t28 = r * r;
    const auto t31 = m * r;
    const auto t34 = t13 * t16;
    const auto t39 = dDdr;
    const auto t52 = dgdtheta;
    const auto t57 = dDdtheta;
    const auto t59 = t14 * (static_cast<scalar_type>(-2) * t13 * t52 + static_cast<scalar_type>(2) * t13 * t20 - t57);
    const auto t62 = t59 / static_cast<scalar_type>(2);
    const auto t66 = t14 * (static_cast<scalar_type>(-2) * t27 + static_cast<scalar_type>(2) * t34 - t39);
    const auto t67 = t66 / static_cast<scalar_type>(2);
    const auto t72 = -t16 * t28 + r + static_cast<scalar_type>(2) * t16 * r * m - m;
    const auto t73 = t25 * t72;
    const auto t76 = st;
    const auto t78 = ct;
    const auto t81 = (t76 * t20 - t78) / t76;
    const auto t83 = t76 * t76;
    const auto t87 = t76 * r;

    christoffel_symbols_type symbols;
    symbols.setZero();
    
    // Non-zero Christoffel symbols
    symbols(0, 0, 1) = -t3 * t7 * t15 * t16;
    symbols(0, 0, 2) = t7 * t12 * t14 * t20;
    
    symbols(0, 1, 0) = t16;
    symbols(0, 2, 0) = t20;
    
    symbols(1, 0, 0) = t16;
    
    symbols(1, 1, 1) = -t25 * t14
        * (static_cast<scalar_type>(2) * t27 * t28 - static_cast<scalar_type>(4) * t27 * t31 
           - static_cast<scalar_type>(2) * t34 * t28 + static_cast<scalar_type>(4) * t34 * t31 
           + t39 * t28 - static_cast<scalar_type>(2) * t39 * r * m
           - static_cast<scalar_type>(2) * t13 * r + static_cast<scalar_type>(2) * t13 * m)
        / static_cast<scalar_type>(2);
    symbols(1, 1, 2) = -t59 * t25 / static_cast<scalar_type>(2);
    
    symbols(1, 2, 1) = -t62;
    symbols(1, 2, 2) = -t67;
    
    symbols(1, 3, 3) = -t73;
    
    symbols(2, 0, 0) = t20;
    
    symbols(2, 1, 1) = -t62;
    symbols(2, 1, 2) = -t67;
    
    symbols(2, 2, 1) = -t3 * t66 / static_cast<scalar_type>(2);
    symbols(2, 2, 2) = -t62;
    
    symbols(2, 3, 3) = -t81;
    
    symbols(3, 1, 3) = -t73;
    
    symbols(3, 2, 3) = -t81;
    
    symbols(3, 3, 1) = t3 * t12 * t14 * t83 * t72;
    symbols(3, 3, 2) = -t15 * t87 * (-t87 * t20 + r * t78 + static_cast<scalar_type>(2) * t76 * m * t20 - static_cast<scalar_type>(2) * m * t78);
    
    return symbols;
  }

  // Mass parameter
  scalar_type mass = static_cast<scalar_type>(1);
  
  // Deformation parameter
  scalar_type q = static_cast<scalar_type>(0.1);
};
}
