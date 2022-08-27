#pragma once

#include <astray/core/metric.hpp>

namespace ast::metrics
{
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class minkowski : public metric<coordinate_system_type::cartesian, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    christoffel_symbols_type symbols;
    symbols.setZero();
    return symbols;
  }
};
}