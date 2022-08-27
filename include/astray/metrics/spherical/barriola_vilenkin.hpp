#pragma once

#include <cmath>

#include <astray/core/metric.hpp>

namespace ast::metrics
{
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class barriola_vilenkin : public metric<coordinate_system_type::spherical, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto t1 = static_cast<scalar_type>(1) / position[1];
    const auto t2 = static_cast<scalar_type>(std::pow(scaling_factor, 2));
    const auto t3 = t2 * position[1];
    const auto t4 = std::sin(position[2]);
    const auto t6 = std::cos(position[2]);
    const auto t7 = static_cast<scalar_type>(1) / t4 * t6;
    const auto t8 = static_cast<scalar_type>(std::pow(t4, 2));

    christoffel_symbols_type symbols;
    symbols.setZero();
    symbols(1, 2, 2) =  t1;
    symbols(1, 3, 3) =  t1;
    symbols(2, 1, 2) =  t1;
    symbols(2, 2, 1) = -t3;
    symbols(2, 3, 3) =  t7;
    symbols(3, 1, 3) =  t1;
    symbols(3, 2, 3) =  t7;
    symbols(3, 3, 1) = -t3 * t8;
    symbols(3, 3, 2) = -t4 * t6;
    return symbols;
  }
  
  scalar_type scaling_factor = static_cast<scalar_type>(1);
};
}