#pragma once

#include <cmath>

#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>

namespace ast::metrics
{
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class friedman_lemaitre_robertson_walker : public metric<coordinate_system_type::spherical, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  enum class curvature_constant : std::int8_t
  {
    negative = -1,
    zero     =  0,
    positive =  1,
  };

  using consts = constants<scalar_type>;
  
  __device__ termination_reason       check_termination  (const vector_type& position, const vector_type& direction) const override
  {
    if (curvature == curvature_constant::negative)
      if (std::abs(static_cast<scalar_type>(1) + static_cast<scalar_type>(0.25) * static_cast<scalar_type>(curvature) * std::pow(position[1], 2)) < consts::epsilon)
        return termination_reason::spacetime_breakdown;
    return termination_reason::none;
  }

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto  k        = static_cast<scalar_type>(curvature);
    const auto  constant = static_cast<scalar_type>(4) * consts::gravitational_constant * mass / consts::three_pi;
    
    scalar_type r ;
    scalar_type dr;
    if      (curvature == curvature_constant::negative)
    {
      r  = static_cast<scalar_type>(0.5 ) * constant * (std::cosh(position[0]) - static_cast<scalar_type>(1));
      dr = static_cast<scalar_type>(0.5 ) * constant *  std::sinh(position[0]);
    }
    else if (curvature == curvature_constant::zero    )
    {
      r  = static_cast<scalar_type>(2.25) * constant * std::pow(position[0], 2);
      dr = static_cast<scalar_type>(4.5 ) * constant *          position[0];
    }
    else if (curvature == curvature_constant::positive)
    {
      r  = static_cast<scalar_type>(0.5 ) * constant * (static_cast<scalar_type>(1) - std::cos(position[0]));
      dr = static_cast<scalar_type>(0.5 ) * constant *                                std::sin(position[0]);
    }

    const auto t2  = static_cast<scalar_type>(1) / r;
    const auto t4  = t2 * dr;
    const auto t5  = std::pow(position[1], 2);
    const auto t6  = k * t5;
    const auto t7  = static_cast<scalar_type>(4) + t6;
    const auto t8  = std::pow(t7, 2);
    const auto t9  = static_cast<scalar_type>(1) / t8;
    const auto t12 = static_cast<scalar_type>(1) / t7;
    const auto t18 = t6 - static_cast<scalar_type>(4);
    const auto t19 = t12 / position[1] * t18;
    const auto t20 = t2 * t5;
    const auto t24 = t12 * position[1];
    const auto t26 = std::sin(position[2]);
    const auto t28 = std::cos(position[2]);
    const auto t29 = static_cast<scalar_type>(1) / t26 * t28;
    const auto t30 = std::pow(t26, 2);

    christoffel_symbols_type symbols;
    symbols.setZero();
    symbols(0, 0, 0) = t4;
    symbols(0, 1, 1) = t4;
    symbols(0, 2, 2) = t4;
    symbols(0, 3, 3) = t4;
    symbols(1, 0, 1) = t4;
    symbols(1, 1, 0) = static_cast<scalar_type>(16) * t4 * t9;
    symbols(1, 1, 1) = static_cast<scalar_type>(-2) * t12 * k * position[1];
    symbols(1, 2, 2) = -t19;
    symbols(1, 3, 3) = -t19;
    symbols(2, 0, 2) = t4;
    symbols(2, 1, 2) = -t19;
    symbols(2, 2, 0) = static_cast<scalar_type>(16) * t20 * t9 * dr;
    symbols(2, 2, 1) = t24 * t18;
    symbols(2, 3, 3) = t29;
    symbols(3, 0, 3) = t4;
    symbols(3, 1, 3) = -t19;
    symbols(3, 2, 3) = t29;
    symbols(3, 3, 0) = static_cast<scalar_type>(16) * t20 * t30 * dr * t9;
    symbols(3, 3, 1) = t24 * t30 * t18;
    symbols(3, 3, 2) = -t26 * t28;
    return symbols;
  }
  
  scalar_type        mass      = static_cast<scalar_type>(1);
  curvature_constant curvature = curvature_constant::positive;
};
}