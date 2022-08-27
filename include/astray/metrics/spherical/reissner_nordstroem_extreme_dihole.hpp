#pragma once

#include <cmath>

#include <astray/core/metric.hpp>

namespace ast::metrics
{
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class reissner_nordstroem_extreme_dihole : public metric<coordinate_system_type::spherical, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto radius1       = std::sqrt(std::pow(position[1], 2) + std::pow(position[2], 2) + std::pow(position[3] - static_cast<scalar_type>(1), 2));
    const auto radius2       = std::sqrt(std::pow(position[1], 2) + std::pow(position[2], 2) + std::pow(position[3] + static_cast<scalar_type>(1), 2));
    const auto u             = static_cast<scalar_type>(1) + mass1 / radius1 + mass2 / radius2;

    const auto radius1_cubed = std::pow(radius1, static_cast<scalar_type>(-3));
    const auto radius2_cubed = std::pow(radius2, static_cast<scalar_type>(-3));
    const auto du            = vector_type {
      static_cast<scalar_type>(0),
      - mass1 * position[1] * radius1_cubed - mass2 * position[1] * radius2_cubed,
      - mass1 * position[2] * radius1_cubed - mass2 * position[2] * radius2_cubed,
      - mass1 * (position[3] - static_cast<scalar_type>(1)) * radius1_cubed 
      - mass2 * (position[3] + static_cast<scalar_type>(1)) * radius2_cubed };

    const auto t2  = static_cast<scalar_type>(std::pow(u , 2));
    const auto t3  = static_cast<scalar_type>(std::pow(t2, 2));
    const auto t5  = static_cast<scalar_type>(1) / t3 / u;
    const auto t12 = static_cast<scalar_type>(1) / u;
    const auto t13 = t12 * du[1];
    const auto t14 = t12 * du[2];
    const auto t15 = t12 * du[3];

    christoffel_symbols_type symbols;
    symbols.setZero();
    symbols(0, 0, 1) = -t5 * du[1];
    symbols(0, 0, 2) = -t5 * du[2];
    symbols(0, 0, 3) = -t5 * du[3];
    symbols(0, 1, 0) = -t13;
    symbols(0, 2, 0) = -t14;
    symbols(0, 3, 0) = -t15;
    symbols(1, 0, 0) = -t13;
    symbols(1, 1, 1) =  t13;
    symbols(1, 1, 2) = -t14;
    symbols(1, 1, 3) = -t15;
    symbols(1, 2, 1) =  t14;
    symbols(1, 2, 2) =  t13;
    symbols(1, 3, 1) =  t15;
    symbols(1, 3, 3) =  t13;
    symbols(2, 0, 0) = -t14;
    symbols(2, 1, 1) =  t14;
    symbols(2, 1, 2) =  t13;
    symbols(2, 2, 1) = -t13;
    symbols(2, 2, 2) =  t14;
    symbols(2, 2, 3) = -t15;
    symbols(2, 3, 2) =  t15;
    symbols(2, 3, 3) =  t14;
    symbols(3, 0, 0) = -t15;
    symbols(3, 1, 1) =  t15;
    symbols(3, 1, 3) =  t13;
    symbols(3, 2, 2) =  t15;
    symbols(3, 2, 3) =  t14;
    symbols(3, 3, 1) = -t13;
    symbols(3, 3, 2) = -t14;
    symbols(3, 3, 3) =  t15;
    return symbols;
  }
  
  scalar_type mass1 = static_cast<scalar_type>(1);
  scalar_type mass2 = static_cast<scalar_type>(1);
};
}