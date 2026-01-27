#pragma once

#include <cmath>

#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>

namespace ast::metrics
{
// Implementation based on Motion4D metric by Thomas Mueller (tauzero7)
// Reference: https://github.com/tauzero7/Motion4D
//
// de Sitter universe in conformal Cartesian coordinates (T, x, y, z).
//
// The metric describes an expanding universe with a cosmological constant.
// In conformal coordinates, the spatial part is conformally flat.
//
// Parameters:
//   hubble - Hubble parameter H
//   p - parameter controlling conformal vs flat (p < 0 gives flat Minkowski)
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class de_sitter_univ_conf : public metric<coordinate_system_type::cartesian, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  __device__ termination_reason check_termination(const vector_type& position, [[maybe_unused]] const vector_type& direction) const override
  {
    // Check for singularity at T = 0 when p >= 0
    if (p >= static_cast<scalar_type>(0) && position[0] <= consts::epsilon)
      return termination_reason::spacetime_breakdown;
    return termination_reason::none;
  }

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto T = position[0];
    auto t1 = static_cast<scalar_type>(1) / T;

    // When p < 0, this reduces to Minkowski spacetime
    if (p < static_cast<scalar_type>(0)) {
      t1 = static_cast<scalar_type>(0);
    }

    christoffel_symbols_type symbols;
    symbols.setZero();
    
    // Non-zero Christoffel symbols (conformal structure)
    symbols(0, 0, 0) = -t1;
    symbols(0, 1, 1) = -t1;
    symbols(0, 2, 2) = -t1;
    symbols(0, 3, 3) = -t1;
    
    symbols(1, 0, 1) = -t1;
    symbols(1, 1, 0) = -t1;
    
    symbols(2, 0, 2) = -t1;
    symbols(2, 2, 0) = -t1;
    
    symbols(3, 0, 3) = -t1;
    symbols(3, 3, 0) = -t1;
    
    return symbols;
  }

  scalar_type hubble = static_cast<scalar_type>(1);
  scalar_type p = static_cast<scalar_type>(1);
};
}
