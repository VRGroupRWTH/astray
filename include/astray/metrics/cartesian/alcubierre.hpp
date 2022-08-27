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
class alcubierre : public metric<coordinate_system_type::cartesian, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto xmvt          = position[1] - velocity * position[0];
    const auto rs            = static_cast<scalar_type>(std::sqrt(
      std::pow(xmvt       , 2) +
      std::pow(position[2], 2) + 
      std::pow(position[3], 2)));

    const auto w1            = std::tanh(thickness * (rs + radius));
    const auto w2            = std::tanh(thickness * (rs - radius));
    const auto w3            = std::tanh(thickness * radius);
    const auto f             = static_cast<scalar_type>(0.5) * (w1 / w3 - w2 / w3);
    
    const auto tanh_positive = std::tanh(thickness * (rs + radius));
    const auto tanh_negative = std::tanh(thickness * (rs - radius));
    const auto factor        = static_cast<scalar_type>(
      (std::pow(tanh_negative, 2) - std::pow(tanh_positive, 2)) * 
      (static_cast<scalar_type>(0.5) * thickness / (rs * std::tanh(thickness * radius))));
    vector_type df {
      -velocity * xmvt * factor,
                  xmvt * factor,
      position[2]      * factor,
      position[3]      * factor};

    const auto t1  = std::pow(velocity, 2);
    const auto t2  = t1 * velocity;
    const auto t3  = f;
    const auto t4  = std::pow(t3, 2);
    const auto t6  = df[1];
    const auto t7  = std::pow(consts::speed_of_light, 2);
    const auto t8  = static_cast<scalar_type>(1) / t7;
    const auto t10 = t2 * t4 * t6 * t8;
    const auto t11 = df[0];
    const auto t14 = t3 * t6;
    const auto t22 = t1 * t3;
    const auto t23 = df[2];
    const auto t25 = df[3];
    const auto t27 = t8 * t1;
    const auto t28 = t27 * t14;
    const auto t29 = velocity * t23;
    const auto t30 = t29                  / static_cast<scalar_type>(2);
    const auto t31 = velocity * t25;
    const auto t32 = t31                  / static_cast<scalar_type>(2);
    const auto t35 = t27 * t3 * t23       / static_cast<scalar_type>(2);
    const auto t38 = (t1 * t4 + t7) * t8;
    const auto t40 = t29 * t38            / static_cast<scalar_type>(2);
    const auto t43 = t27 * t3 * t25       / static_cast<scalar_type>(2);
    const auto t45 = t31 * t38            / static_cast<scalar_type>(2);
    const auto t46 = velocity * t8;
    const auto t49 = t46 * t23            / static_cast<scalar_type>(2);
    const auto t51 = t46 * t25            / static_cast<scalar_type>(2);

    christoffel_symbols_type symbols;
    symbols.setZero();
    symbols(0, 0, 0) =  t10;
    symbols(0, 0, 1) =  velocity * (-t7 * t11 - t7 * velocity * t14 + t2 * t4 * t3 * t6) * t8;
    symbols(0, 0, 2) = -t22 * t23;
    symbols(0, 0, 3) = -t22 * t25;
    symbols(0, 1, 0) = -t28;
    symbols(0, 1, 1) = -t10;
    symbols(0, 1, 2) =  t30;
    symbols(0, 1, 3) =  t32;
    symbols(0, 2, 0) = -t35;
    symbols(0, 2, 1) = -t40;
    symbols(0, 3, 0) = -t43;
    symbols(0, 3, 1) = -t45;
    symbols(1, 0, 0) = -t28;
    symbols(1, 0, 1) = -t10;
    symbols(1, 0, 2) =  t30;
    symbols(1, 0, 3) =  t32;
    symbols(1, 1, 0) =  t46 * t6;
    symbols(1, 1, 1) =  t28;
    symbols(1, 2, 0) =  t49;
    symbols(1, 2, 1) =  t35;
    symbols(1, 3, 0) =  t51;
    symbols(1, 3, 1) =  t43;
    symbols(2, 0, 0) = -t35;
    symbols(2, 0, 1) = -t40;
    symbols(2, 1, 0) =  t49;
    symbols(2, 1, 1) =  t35;
    symbols(3, 0, 0) = -t43;
    symbols(3, 0, 1) = -t45;
    symbols(3, 1, 0) =  t51;
    symbols(3, 1, 1) =  t43;
    return symbols;
  }

  scalar_type velocity  = static_cast<scalar_type>(0);
  scalar_type radius    = static_cast<scalar_type>(1);
  scalar_type thickness = static_cast<scalar_type>(1);
};
}