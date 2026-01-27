#pragma once

#include <cmath>

#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>

namespace ast::metrics
{
// Implementation based on Motion4D metric by Thomas Mueller (tauzero7)
// Reference: https://github.com/tauzero7/Motion4D
//
// Incoming radiation Vaidya metric in spherical coordinates (v,r,theta,phi).
//
// The line element is given by:
//   ds^2 = 2 dv dr - (1 - 2m(v)/r) dv^2 + r^2 (dθ^2 + sin²θ dφ^2)
//
// This represents a black hole with time-varying mass m(v) due to incoming radiation.
// The mass function is m(v) = k * (vl - v)^(1/3) * tanh²(σ(vl - v)) for v < vl, and 0 for v >= vl.
//
// Reference: Piesnack and Kassner, https://arxiv.org/pdf/2103.08340v2.pdf
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class vaidya_incoming_radiation : public metric<coordinate_system_type::spherical, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  __device__ termination_reason check_termination(const vector_type& position, [[maybe_unused]] const vector_type& direction) const override
  {
    const auto v = position[0];
    const auto r = position[1];
    scalar_type m;
    calc_mass_function(v, m);
    
    // Check for singularity - using same condition as Motion4D
    if (r * r <= (static_cast<scalar_type>(1) + consts::epsilon) * static_cast<scalar_type>(4) * m * m)
      return termination_reason::spacetime_breakdown;
    return termination_reason::none;
  }

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto r     = position[1];
    const auto theta = position[2];
    
    scalar_type m, dmdv;
    calc_mass_function(position[0], m, dmdv);
    
    const auto t1  = m;
    const auto t2  = r * r;
    const auto t4  = t1 / t2;
    const auto t5  = dmdv;
    const auto t8  = t1 * t1;
    const auto t14 = static_cast<scalar_type>(1) / r;
    const auto t16 = -r + static_cast<scalar_type>(2) * t1;
    const auto t17 = std::sin(theta);
    const auto t19 = std::cos(theta);
    const auto t20 = static_cast<scalar_type>(1) / t17 * t19;
    const auto t21 = t17 * t17;

    christoffel_symbols_type symbols;
    symbols.setZero();
    
    // Non-zero Christoffel symbols
    symbols(0, 0, 0) = t4;
    symbols(0, 0, 1) = -(-t5 * t2 - t1 * r + static_cast<scalar_type>(2) * t8) / t2 / r;
    
    symbols(0, 1, 1) = -t4;
    
    symbols(1, 0, 1) = -t4;
    symbols(1, 2, 2) = t14;
    symbols(1, 3, 3) = t14;
    
    symbols(2, 1, 2) = t14;
    symbols(2, 2, 0) = -r;
    symbols(2, 2, 1) = t16;
    symbols(2, 3, 3) = t20;
    
    symbols(3, 1, 3) = t14;
    symbols(3, 2, 3) = t20;
    symbols(3, 3, 0) = -r * t21;
    symbols(3, 3, 1) = t16 * t21;
    symbols(3, 3, 2) = -t17 * t19;
    
    return symbols;
  }

private:
  __device__ void calc_mass_function(scalar_type v, scalar_type& m) const
  {
    if (v >= vl) {
      m = static_cast<scalar_type>(0);
      return;
    }
    
    const auto th = std::tanh(sigma * (vl - v));
    m = k * std::pow(vl - v, static_cast<scalar_type>(1) / static_cast<scalar_type>(3)) * th * th;
  }

  __device__ void calc_mass_function(scalar_type v, scalar_type& m, scalar_type& dmdv) const
  {
    if (v >= vl) {
      m = dmdv = static_cast<scalar_type>(0);
      return;
    }
    
    const auto th = std::tanh(sigma * (vl - v));
    const auto ch = std::cosh(sigma * (vl - v));
    
    m = k * std::pow(vl - v, static_cast<scalar_type>(1) / static_cast<scalar_type>(3)) * th * th;
    
    dmdv = -k / static_cast<scalar_type>(3) * std::pow(vl - v, -static_cast<scalar_type>(2) / static_cast<scalar_type>(3)) * th * th;
    dmdv -= k * std::pow(vl - v, static_cast<scalar_type>(1) / static_cast<scalar_type>(3)) * static_cast<scalar_type>(2) * sigma * th / (ch * ch);
  }

public:
  scalar_type k     = static_cast<scalar_type>(1);      // Constant factor
  scalar_type vl    = static_cast<scalar_type>(0);      // v limit
  scalar_type sigma = static_cast<scalar_type>(0.5);    // Sigma parameter
};
}
