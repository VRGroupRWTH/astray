#pragma once

#include <cmath>

#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>
#include <astray/math/ipow.hpp>

namespace ast::metrics
{
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class kottler : public metric<coordinate_system_type::spherical, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;
  
  __device__ termination_reason       check_termination  (const vector_type& position, const vector_type& direction) const override
  {
    const auto rs = consts::schwarzschild_radius(mass);
    const auto r = position[1];
    if (r < static_cast<scalar_type>(0))
      return termination_reason::spacetime_breakdown;
    const auto r_sq = ipow<2>(r);
    if (std::abs(static_cast<scalar_type>(1) - rs / r 
      - consts::cosmological_constant / static_cast<scalar_type>(3) * r_sq) <= consts::epsilon)
      return termination_reason::spacetime_breakdown;
    return termination_reason::none;
  }

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto rs  = consts::schwarzschild_radius(mass);
    const auto r   = position[1];
    const auto r_sq = ipow<2>(r);

    const auto t2  = static_cast<scalar_type>(3) * rs;
    const auto t4  = r_sq * r;
    const auto t5  = consts::cosmological_constant * t4;
    const auto t6  = -static_cast<scalar_type>(3) * r + t2 + t5;
    const auto t11 = -t2 + static_cast<scalar_type>(2) * t5;
    const auto t15 = static_cast<scalar_type>(1) / r;
    const auto t19 = t15 / t6 * t11 / static_cast<scalar_type>(2);
    const auto t22 = std::sin(position[2]);
    const auto t24 = std::cos(position[2]);
    const auto t25 = static_cast<scalar_type>(1) / t22 * t24;
    const auto t22_sq = ipow<2>(t22);

    christoffel_symbols_type symbols;
    symbols.setZero();
    symbols(0, 0, 1) =  t6 / t4 * consts::speed_of_light_squared * t11 / static_cast<scalar_type>(18);
    symbols(0, 1, 0) =  t19;
    symbols(1, 0, 0) =  t19;
    symbols(1, 1, 1) = -t19;
    symbols(1, 2, 2) =  t15;
    symbols(1, 3, 3) =  t15;
    symbols(2, 1, 2) =  t15;
    symbols(2, 2, 1) = -r + rs + t5 / static_cast<scalar_type>(3);
    symbols(2, 3, 3) =  t25;
    symbols(3, 1, 3) =  t15;
    symbols(3, 2, 3) =  t25;
    symbols(3, 3, 1) =  t6  * t22_sq / static_cast<scalar_type>(3);
    symbols(3, 3, 2) = -t22 * t24;
    return symbols;
  }

  scalar_type mass = static_cast<scalar_type>(1);
};
}