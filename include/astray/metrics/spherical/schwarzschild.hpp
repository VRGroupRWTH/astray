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
class schwarzschild : public metric<coordinate_system_type::spherical, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  __device__ termination_reason       check_termination  (const vector_type& position, const vector_type& direction) const override
  {
    const auto rs = consts::schwarzschild_radius(mass);
    if (position[1] < static_cast<scalar_type>(0) || 
        static_cast<scalar_type>(std::pow(position[1], 2)) <= (static_cast<scalar_type>(1) + consts::epsilon) * static_cast<scalar_type>(std::pow(rs, 2)))
      return termination_reason::spacetime_breakdown;
    return termination_reason::none;
  }

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto rs    = consts::schwarzschild_radius(mass);
    const auto r     = position[1];
    const auto theta = position[2];
    const auto t1    = r - rs;
    const auto t2    = static_cast<scalar_type>(std::pow(r, 2));
    const auto t6    = consts::speed_of_light_squared;
    const auto t10   = static_cast<scalar_type>(1) / r;
    const auto t14   = t10 / t1 * rs * static_cast<scalar_type>(0.5);
    const auto t15   = std::sin(theta);
    const auto t17   = std::cos(theta);
    const auto t18   = static_cast<scalar_type>(1) / t15 * t17;
    const auto t19   = static_cast<scalar_type>(std::pow(t15, 2));

    christoffel_symbols_type symbols;
    symbols.setZero();
    symbols(0, 0, 1) = t1 / t2 / r * t6 * rs / static_cast<scalar_type>(2);
    symbols(0, 1, 0) =  t14;
    symbols(1, 0, 0) =  t14;
    symbols(1, 1, 1) = -t14;
    symbols(1, 2, 2) =  t10;
    symbols(1, 3, 3) =  t10;
    symbols(2, 1, 2) =  t10;
    symbols(2, 2, 1) = -t1;
    symbols(2, 3, 3) =  t18;
    symbols(3, 1, 3) =  t10;
    symbols(3, 2, 3) =  t18;
    symbols(3, 3, 1) = -t1  * t19;
    symbols(3, 3, 2) = -t15 * t17;
    return symbols;
  }
  
  scalar_type mass = static_cast<scalar_type>(1);
};
}