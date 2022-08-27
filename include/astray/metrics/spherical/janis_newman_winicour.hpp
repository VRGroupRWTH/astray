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
class janis_newman_winicour : public metric<coordinate_system_type::spherical, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  __device__ termination_reason       check_termination  (const vector_type& position, const vector_type& direction) const override
  {
    const auto rs = consts::schwarzschild_radius(mass);
    if (position[1] < static_cast<scalar_type>(0) || 
        position[1] < rs * static_cast<scalar_type>(std::pow(static_cast<scalar_type>(0.5) * (static_cast<scalar_type>(1) + gamma) / gamma, 2)) ||
        static_cast<scalar_type>(std::pow(gamma, 2) * std::pow(position[1], 2)) <= (static_cast<scalar_type>(1) + consts::epsilon) * static_cast<scalar_type>(std::pow(rs, 2)))
      return termination_reason::spacetime_breakdown;
    return termination_reason::none;
  }

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto rs  = consts::schwarzschild_radius(mass);

    const auto t1  = gamma * position[1];
    const auto t2  = t1 - rs;
    const auto t3  = static_cast<scalar_type>(1) / gamma;
    const auto t5  = static_cast<scalar_type>(1) / position[1];
    const auto t7  = static_cast<scalar_type>(std::pow(t2 * t3 * t5, gamma));
    const auto t8  = static_cast<scalar_type>(std::pow(t7, 2));
    const auto t13 = static_cast<scalar_type>(1) / t2;
    const auto t17 = rs * gamma;
    const auto t18 = t5 * t13;
    const auto t20 = t17 * t18 / static_cast<scalar_type>(2);
    const auto t22 = static_cast<scalar_type>(2) * t1 - t17 - rs;
    const auto t24 = t18 * t22 / static_cast<scalar_type>(2);
    const auto t27 = std::sin(position[2]);
    const auto t29 = std::cos(position[2]);
    const auto t30 = static_cast<scalar_type>(1) / t27 * t29;
    const auto t31 = static_cast<scalar_type>(std::pow(t27, 2));

    christoffel_symbols_type symbols;
    symbols.setZero();
    symbols(0, 0, 1) =  t8 * rs * gamma * t5 * consts::speed_of_light_squared * t13 / static_cast<scalar_type>(2);
    symbols(0, 1, 0) =  t20;
    symbols(1, 0, 0) =  t20;
    symbols(1, 1, 1) = -t20;
    symbols(1, 2, 2) =  t24;
    symbols(1, 3, 3) =  t24;
    symbols(2, 1, 2) =  t24;
    symbols(2, 2, 1) = -t22 * t3 / static_cast<scalar_type>(2);
    symbols(2, 3, 3) =  t30;
    symbols(3, 1, 3) =  t24;
    symbols(3, 2, 3) =  t30;
    symbols(3, 3, 1) = -t22 * t31 * t3 / static_cast<scalar_type>(2);
    symbols(3, 3, 2) = -t27 * t29;
    return symbols;
  }
  
  scalar_type mass  = static_cast<scalar_type>(1);
  scalar_type gamma = static_cast<scalar_type>(1);
};
}