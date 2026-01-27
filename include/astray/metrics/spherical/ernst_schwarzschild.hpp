#pragma once

#include <cmath>

#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>

namespace ast::metrics
{
// Implementation based on Motion4D metric by Thomas Mueller (tauzero7)
// Reference: https://github.com/tauzero7/Motion4D
//
// Ernst-Schwarzschild metric in spherical coordinates (t,r,theta,phi).
//
// The line element is given by:
//   ds² = Λ²[(1-2m/r) dt² - dr²/(1-2m/r) - r²dθ²] - (r²sin²θ/Λ²)dφ²
// where Λ = 1 + B²r²sin²θ.
//
// This is the Ernst form of the Schwarzschild metric in a magnetic universe.
//
// References:
// - Frederick J. Ernst, "Black holes in a magnetic universe," J. Math. Phys. 17, 54-56 (1976).
// - R.A. Konoplya, "Magnetised black hole as a gravitational lens," Phys. Lett. B 644, 219-223 (2007).
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class ernst_schwarzschild : public metric<coordinate_system_type::spherical, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  __device__ termination_reason check_termination(const vector_type& position, [[maybe_unused]] const vector_type& direction) const override
  {
    const auto r = position[1];
    const auto rs = static_cast<scalar_type>(2) * mass;
    
    // Check for singularity at r = 0
    if (r <= static_cast<scalar_type>(0))
      return termination_reason::spacetime_breakdown;
    
    // Check for event horizon at r = rs
    if (r <= rs)
      return termination_reason::event_horizon;
    
    return termination_reason::none;
  }

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto r = position[1];
    const auto theta = position[2];
    const auto m = mass;
    const auto B = magnetic_field;

    const auto t1 = B * B;
    const auto t2 = r * r;
    const auto t3 = t2 * r;
    const auto t5 = std::sin(theta);
    const auto t6 = t5 * t5;
    const auto t8 = static_cast<scalar_type>(2) * t1 * t3 * t6;
    const auto t9 = t1 * t2;
    const auto t11 = t9 * t6 * m;
    const auto t13 = t8 - static_cast<scalar_type>(3) * t11 + m;
    const auto t15 = r - static_cast<scalar_type>(2) * m;
    const auto t18 = t9 * t6;
    const auto t19 = static_cast<scalar_type>(1) + t18;
    const auto t20 = static_cast<scalar_type>(1) / t19;
    const auto t24 = static_cast<scalar_type>(1) / r;
    const auto t26 = std::cos(theta);
    const auto t28 = t5 * t26 * t20;
    const auto t31 = static_cast<scalar_type>(1) / t15;
    const auto t33 = t24 * t20;
    const auto t34 = t13 * t31 * t33;
    const auto t36 = static_cast<scalar_type>(2) * t9 * t28;
    const auto t48 = static_cast<scalar_type>(3) * t18 + static_cast<scalar_type>(1);
    const auto t50 = t48 * t24 * t20;
    const auto t51 = -static_cast<scalar_type>(1) + t18;
    const auto t53 = t51 * t24 * t20;
    const auto t56 = t51 * t26;
    const auto t59 = t56 * t20 / t5;
    const auto t61 = t1 * t1;
    const auto t62 = t2 * t2;
    const auto t64 = t6 * t6;
    const auto t67 = static_cast<scalar_type>(1) / (static_cast<scalar_type>(1) + static_cast<scalar_type>(2) * t18 + t61 * t62 * t64);
    const auto t70 = t19 * t19;
    const auto t72 = static_cast<scalar_type>(1) / t70 / t19;

    christoffel_symbols_type symbols;
    symbols.setZero();
    
    // Non-zero Christoffel symbols
    symbols(0, 0, 1) = t13 * t15 / t3 * t20;
    symbols(0, 0, 2) = static_cast<scalar_type>(2) * t15 * t1 * t24 * t28;
    symbols(0, 1, 0) = t34;
    symbols(0, 2, 0) = t36;
    
    symbols(1, 0, 0) = t34;
    symbols(1, 1, 1) = (t8 - static_cast<scalar_type>(5) * t11 - m) * t31 * t33;
    symbols(1, 1, 2) = -static_cast<scalar_type>(2) * t1 * r * t5 * t26 * t31 * t20;
    symbols(1, 2, 1) = t36;
    symbols(1, 2, 2) = t50;
    symbols(1, 3, 3) = -t53;
    
    symbols(2, 0, 0) = t36;
    symbols(2, 1, 1) = t36;
    symbols(2, 1, 2) = t50;
    symbols(2, 2, 1) = -t48 * t15 * t20;
    symbols(2, 2, 2) = t36;
    symbols(2, 3, 3) = -t59;
    
    symbols(3, 1, 3) = -t53;
    symbols(3, 2, 3) = -t59;
    symbols(3, 3, 1) = t15 * t67 * t6 * t51 * t72;
    symbols(3, 3, 2) = t67 * t5 * t56 * t72;
    
    return symbols;
  }

  scalar_type mass = static_cast<scalar_type>(1);
  scalar_type magnetic_field = static_cast<scalar_type>(0.1);
};
}
