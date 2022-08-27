#pragma once

#include <cmath>

#include <astray/core/metric.hpp>

namespace ast::metrics
{
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class morris_thorne : public metric<coordinate_system_type::spherical, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto t1  = static_cast<scalar_type>(std::pow(position[1]  , 2));
    const auto t2  = static_cast<scalar_type>(std::pow(throat_radius, 2));
    const auto t5  = static_cast<scalar_type>(1) / (t1 + t2) * position[1];
    const auto t6  = std::sin(position[2]);
    const auto t8  = std::cos(position[2]);
    const auto t9  = static_cast<scalar_type>(1) / t6 * t8;
    const auto t10 = static_cast<scalar_type>(std::pow(t6, 2));

    christoffel_symbols_type symbols;
    symbols.setZero();
    symbols(1, 2, 2) = t5;
    symbols(1, 3, 3) = t5;
    symbols(2, 1, 2) = t5;
    symbols(2, 2, 1) = -position[1];
    symbols(2, 3, 3) = t9;
    symbols(3, 1, 3) = t5;
    symbols(3, 2, 3) = t9;
    symbols(3, 3, 1) = -position[1] * t10;
    symbols(3, 3, 2) = -t6 * t8;
    return symbols;
  }
  
  scalar_type throat_radius = static_cast<scalar_type>(1);
};
}