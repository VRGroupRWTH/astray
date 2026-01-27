#pragma once

#include <cmath>

#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>

namespace ast::metrics
{
// Implementation based on Motion4D metric by Thomas Mueller (tauzero7)
// Reference: https://github.com/tauzero7/Motion4D
//
// Sultana-Dyer cosmological black hole in spherical coordinates (t,r,theta,phi).
//
// The line element is given by:
//   ds² = t⁴[(1-2m/r)dt² - 4m/r dt dr - (1+2m/r)dr² - r²(dθ² + sin²θ dφ²)]
//
// This describes a black hole in an Einstein-de Sitter universe.
//
// References:
// - J. Sultana and C.C. Dyer, "Cosmological black holes: A black hole in the Einstein-de Sitter universe",
//   Gen. Relativ. Gravit. 37, 1349-1370 (2005).
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class sultana_dyer : public metric<coordinate_system_type::spherical, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  __device__ termination_reason check_termination(const vector_type& position, [[maybe_unused]] const vector_type& direction) const override
  {
    const auto r = position[1];
    
    // Check for singularity at r = 0
    if (r <= static_cast<scalar_type>(0))
      return termination_reason::spacetime_breakdown;
    
    return termination_reason::none;
  }

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto t = position[0];
    const auto r = position[1];
    const auto theta = position[2];
    const auto m = mass;

    const auto t1 = r * r;
    const auto t2 = t1 * r;
    const auto t3 = m * m;
    const auto t5 = static_cast<scalar_type>(4) * r * t3;
    const auto t6 = t3 * t;
    const auto t8 = static_cast<scalar_type>(1) / t;
    const auto t10 = static_cast<scalar_type>(1) / t2;
    const auto t13 = static_cast<scalar_type>(2) * m;
    const auto t19 = (static_cast<scalar_type>(4) * r + t) * t8 * t10;
    const auto t23 = (r + t13) * m * t19;
    const auto t27 = static_cast<scalar_type>(2) * (t6 - t2 + t5) * t8 * t10;
    const auto t28 = static_cast<scalar_type>(2) * t8;
    const auto t29 = t * r;
    const auto t37 = t * m;
    const auto t38 = static_cast<scalar_type>(2) * t37;
    const auto t40 = r * m;
    const auto t46 = static_cast<scalar_type>(1) / r;
    const auto t48 = -t1 - static_cast<scalar_type>(2) * t40 + t37;
    const auto t52 = static_cast<scalar_type>(4) * t40 + t29 - t38;
    const auto t54 = std::sin(theta);
    const auto t56 = std::cos(theta);
    const auto t57 = static_cast<scalar_type>(1) / t54 * t56;
    const auto t58 = t54 * t54;

    christoffel_symbols_type symbols;
    symbols.setZero();
    
    // Non-zero Christoffel symbols
    symbols(0, 0, 0) = static_cast<scalar_type>(2) * (t2 + t5 + t6) * t8 * t10;
    symbols(0, 0, 1) = m * (r - t13) * t19;
    symbols(0, 1, 0) = t23;
    symbols(0, 1, 1) = -t27;
    symbols(0, 2, 2) = t28;
    symbols(0, 3, 3) = t28;
    
    symbols(1, 0, 0) = t23;
    symbols(1, 0, 1) = -t27;
    symbols(1, 1, 0) = static_cast<scalar_type>(2) * (t29 * m + t2 + static_cast<scalar_type>(4) * t1 * m + t6 + t5) * t8 * t10;
    symbols(1, 1, 1) = -m * (t38 + static_cast<scalar_type>(4) * t1 + static_cast<scalar_type>(8) * t40 + t29) * t8 * t10;
    symbols(1, 2, 2) = t46;
    symbols(1, 3, 3) = t46;
    
    symbols(2, 0, 2) = t28;
    symbols(2, 1, 2) = t46;
    symbols(2, 2, 0) = -static_cast<scalar_type>(2) * t48 * t8;
    symbols(2, 2, 1) = -t52 * t8;
    symbols(2, 3, 3) = t57;
    
    symbols(3, 0, 3) = t28;
    symbols(3, 1, 3) = t46;
    symbols(3, 2, 3) = t57;
    symbols(3, 3, 0) = -static_cast<scalar_type>(2) * t58 * t48 * t8;
    symbols(3, 3, 1) = -t58 * t52 * t8;
    symbols(3, 3, 2) = -t54 * t56;
    
    return symbols;
  }

  scalar_type mass = static_cast<scalar_type>(1);
};
}
