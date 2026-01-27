#pragma once

#include <cmath>

#include <astray/core/metric.hpp>

namespace ast::metrics
{
// Implementation based on Motion4D metric by Thomas Mueller (tauzero7)
// Reference: https://github.com/tauzero7/Motion4D
//
// Minkowski metric in conformal spherical coordinates (ψ, ξ, θ, φ).
//
// The line element is given by:
//   ds^2 = -dψ^2 + dξ^2 + sin²ξ (dθ^2 + sin²θ dφ^2)
//
// These coordinates represent flat Minkowski spacetime in a conformally
// compactified form, which is useful for studying causal structure.
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class minkowski_conformal : public metric<coordinate_system_type::spherical, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto xi    = position[1];
    const auto theta = position[2];
    
    const auto t1  = std::sin(xi);
    const auto t3  = std::cos(xi);
    const auto t4  = static_cast<scalar_type>(1) / t1 * t3;
    const auto t6  = std::sin(theta);
    const auto t8  = std::cos(theta);
    const auto t9  = static_cast<scalar_type>(1) / t6 * t8;
    const auto t10 = t6 * t6;

    christoffel_symbols_type symbols;
    symbols.setZero();
    
    // Non-zero Christoffel symbols
    symbols(1, 2, 2) = t4;
    symbols(1, 3, 3) = t4;
    
    symbols(2, 1, 2) = t4;
    symbols(2, 2, 1) = -t1 * t3;
    symbols(2, 3, 3) = t9;
    
    symbols(3, 1, 3) = t4;
    symbols(3, 2, 3) = t9;
    symbols(3, 3, 1) = -t1 * t10 * t3;
    symbols(3, 3, 2) = -t6 * t8;
    
    return symbols;
  }
};
}
