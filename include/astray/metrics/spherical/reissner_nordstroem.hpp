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
class reissner_nordstroem : public metric<coordinate_system_type::spherical, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  __device__ termination_reason       check_termination  (const vector_type& position, const vector_type& direction) const override
  {
    if (position[1] <= static_cast<scalar_type>(0))
      return termination_reason::numeric_error;
    
    const auto half_schwarzschild_radius = static_cast<scalar_type>(0.5) * consts::schwarzschild_radius(mass);
    const auto disk_radius               = static_cast<scalar_type>(1) - static_cast<scalar_type>(4) * 
      consts::characteristic_length_scale(charge) / std::pow(consts::schwarzschild_radius(mass), 2);
    const auto sqrt_disk_radius          = std::sqrt(disk_radius);
    const auto outer_horizon             = disk_radius < 0 ? static_cast<scalar_type>(0) : half_schwarzschild_radius * (static_cast<scalar_type>(1) + sqrt_disk_radius);
    const auto inner_horizon             = disk_radius < 0 ? static_cast<scalar_type>(0) : half_schwarzschild_radius * (static_cast<scalar_type>(1) - sqrt_disk_radius);

    if (disk_radius > static_cast<scalar_type>(0) && position[1] >= inner_horizon && position[1] <= outer_horizon)
      return termination_reason::spacetime_breakdown;

    return termination_reason::none;
  }

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto t1  = std::pow(position[1], 2);
    const auto t2  = consts::schwarzschild_radius(mass) * position[1];
    const auto t3  = std::pow(charge, 2);
    const auto t4  = consts::characteristic_length_scale(charge);
    const auto t5  = t1 - t2 + t4;
    const auto t6  = std::pow(t1, 2);
    const auto t10 = consts::speed_of_light_squared;
    const auto t12 = t2 - static_cast<scalar_type>(2) * t4;
    const auto t16 = static_cast<scalar_type>(1) / position[1];
    const auto t20 = t16 / t5 * t12 / static_cast<scalar_type>(2);
    const auto t21 = t5 * t16;
    const auto t22 = std::sin(position[2]);
    const auto t24 = std::cos(position[2]);
    const auto t25 = static_cast<scalar_type>(1) / t22 * t24;
    const auto t26 = std::pow(t22, 2);

    christoffel_symbols_type symbols;
    symbols.setZero();
    symbols(0, 0, 1) =  t5 / t6 / position[1] * t10 * t12 / static_cast<scalar_type>(2);
    symbols(0, 1, 0) =  t20;
    symbols(1, 0, 0) =  t20;
    symbols(1, 1, 1) = -t20;
    symbols(1, 2, 2) =  t16;
    symbols(1, 3, 3) =  t16;
    symbols(2, 1, 2) =  t16;
    symbols(2, 2, 1) = -t21;
    symbols(2, 3, 3) =  t25;
    symbols(3, 1, 3) =  t16;
    symbols(3, 2, 3) =  t25;
    symbols(3, 3, 1) = -t21 * t26;
    symbols(3, 3, 2) = -t22 * t24;
    return symbols;
  }
  
  scalar_type mass   = static_cast<scalar_type>(1);
  scalar_type charge = static_cast<scalar_type>(1);
};
}