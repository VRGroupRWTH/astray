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
class de_sitter : public metric<coordinate_system_type::cartesian, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto t1 = consts::speed_of_light_squared;
    const auto t4 = std::exp(hubble_parameter * position[0]);
    const auto t5 = std::pow(t4, 2);
    const auto t7 = static_cast<scalar_type>(1) / t1 * t5 * hubble_parameter;

    christoffel_symbols_type symbols;
    symbols.setZero();
    symbols(0, 1, 1) = hubble_parameter;
    symbols(0, 2, 2) = hubble_parameter;
    symbols(0, 3, 3) = hubble_parameter;
    symbols(1, 0, 1) = hubble_parameter;
    symbols(1, 1, 0) = t7;
    symbols(2, 0, 2) = hubble_parameter;
    symbols(2, 2, 0) = t7;
    symbols(3, 0, 3) = hubble_parameter;
    symbols(3, 3, 0) = t7;
    return symbols;
  }

  scalar_type hubble_parameter = std::sqrt(consts::cosmological_constant * consts::speed_of_light_squared / static_cast<scalar_type>(3));
};
}