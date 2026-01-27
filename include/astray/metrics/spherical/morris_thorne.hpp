#pragma once

#include <cmath>

#include <astray/core/metric.hpp>
#include <astray/math/ipow.hpp>

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
    const auto r   = position[1];
    const auto r_sq = ipow<2>(r);
    const auto b_sq = ipow<2>(throat_radius);
    const auto t5  = static_cast<scalar_type>(1) / (r_sq + b_sq) * r;
    const auto t6  = std::sin(position[2]);
    const auto t8  = std::cos(position[2]);
    const auto t9  = static_cast<scalar_type>(1) / t6 * t8;
    const auto t6_sq = ipow<2>(t6);

    christoffel_symbols_type symbols;
    symbols.setZero();
    symbols(1, 2, 2) = t5;
    symbols(1, 3, 3) = t5;
    symbols(2, 1, 2) = t5;
    symbols(2, 2, 1) = -r;
    symbols(2, 3, 3) = t9;
    symbols(3, 1, 3) = t5;
    symbols(3, 2, 3) = t9;
    symbols(3, 3, 1) = -r * t6_sq;
    symbols(3, 3, 2) = -t6 * t8;
    return symbols;
  }
  
  scalar_type throat_radius = static_cast<scalar_type>(1);
};
}