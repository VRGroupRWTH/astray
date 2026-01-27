#pragma once

#include <cmath>

#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>

namespace ast::metrics
{
// Implementation based on Motion4D metric by Thomas Mueller (tauzero7)
// Reference: https://github.com/tauzero7/Motion4D
//
// Ernst metric in spherical Schwarzschild-like coordinates (t,r,theta,phi).
//
// The line element is given by:
//   ds² = Λ²[-(1-2m/r) dt² + dr²/(1-2m/r) + r²dθ²] + (r²sin²θ/Λ²)dφ²
// where Λ = 1 + B²r²sin²θ.
//
// This describes a black hole in a magnetic universe.
//
// References:
// - Frederick J. Ernst, "Black holes in a magnetic universe," J. Math. Phys. 17, 54-56 (1976).
// - R.A. Konoplya, "Magnetised black hole as a gravitational lens," Phys. Lett. B 644, 219-223 (2007).
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class ernst : public metric<coordinate_system_type::spherical, scalar_type, vector_type, christoffel_symbols_type>
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
    const auto c = consts::speed_of_light;
    const auto rs = static_cast<scalar_type>(2) * mass;
    const auto B = magnetic_field;

    const auto t1 = B * B;
    const auto t2 = r * r;
    const auto t3 = t2 * r;
    const auto t5 = std::sin(theta);
    const auto t6 = t5 * t5;
    const auto t8 = static_cast<scalar_type>(4) * t1 * t3 * t6;
    const auto t11 = t1 * t6 * rs * t2;
    const auto t13 = t8 - static_cast<scalar_type>(3) * t11 + rs;
    const auto t14 = r - rs;
    const auto t16 = c * c;
    const auto t17 = t1 * t2;
    const auto t18 = t17 * t6;
    const auto t19 = static_cast<scalar_type>(1) + t18;
    const auto t20 = static_cast<scalar_type>(1) / t19;
    const auto t27 = static_cast<scalar_type>(1) / r;
    const auto t30 = std::cos(theta);
    const auto t35 = static_cast<scalar_type>(1) / t14;
    const auto t37 = t27 * t20;
    const auto t39 = t13 * t35 * t37 / static_cast<scalar_type>(2);
    const auto t40 = t5 * t30;
    const auto t43 = static_cast<scalar_type>(2) * t17 * t40 * t20;
    const auto t55 = static_cast<scalar_type>(3) * t18 + static_cast<scalar_type>(1);
    const auto t57 = t55 * t27 * t20;
    const auto t58 = -static_cast<scalar_type>(1) + t18;
    const auto t60 = t58 * t27 * t20;
    const auto t63 = t58 * t30;
    const auto t66 = t63 * t20 / t5;
    const auto t68 = t1 * t1;
    const auto t69 = t2 * t2;
    const auto t71 = t6 * t6;
    const auto t74 = static_cast<scalar_type>(1) / (static_cast<scalar_type>(1) + static_cast<scalar_type>(2) * t18 + t68 * t69 * t71);
    const auto t77 = t19 * t19;
    const auto t79 = static_cast<scalar_type>(1) / t77 / t19;

    christoffel_symbols_type symbols;
    symbols.setZero();
    
    // Non-zero Christoffel symbols
    symbols(0, 0, 1) = t13 * t14 * t16 * t20 / t3 / static_cast<scalar_type>(2);
    symbols(0, 0, 2) = static_cast<scalar_type>(2) * t14 * t1 * t27 * t5 * t16 * t20 * t30;
    symbols(0, 1, 0) = t39;
    symbols(0, 2, 0) = t43;
    
    symbols(1, 0, 0) = t39;
    symbols(1, 1, 1) = (t8 - static_cast<scalar_type>(5) * t11 - rs) * t35 * t37 / static_cast<scalar_type>(2);
    symbols(1, 1, 2) = -static_cast<scalar_type>(2) * t40 * r * t1 * t35 * t20;
    symbols(1, 2, 1) = t43;
    symbols(1, 2, 2) = t57;
    symbols(1, 3, 3) = -t60;
    
    symbols(2, 0, 0) = t43;
    symbols(2, 1, 1) = t43;
    symbols(2, 1, 2) = t57;
    symbols(2, 2, 1) = -t55 * t14 * t20;
    symbols(2, 2, 2) = t43;
    symbols(2, 3, 3) = -t66;
    
    symbols(3, 1, 3) = -t60;
    symbols(3, 2, 3) = -t66;
    symbols(3, 3, 1) = t14 * t74 * t6 * t58 * t79;
    symbols(3, 3, 2) = t74 * t5 * t63 * t79;
    
    return symbols;
  }

  scalar_type mass = static_cast<scalar_type>(1);
  scalar_type magnetic_field = static_cast<scalar_type>(0.1);
};
}
