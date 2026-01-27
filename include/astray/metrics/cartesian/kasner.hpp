#pragma once

#include <cmath>

#include <astray/core/metric.hpp>

namespace ast::metrics
{
// Implementation based on Motion4D metric by Thomas Mueller (tauzero7)
// Reference: https://github.com/tauzero7/Motion4D
// 
// Kasner metric in cartesian coordinates (t,x,y,z).
//
// The line element is given by:
//   ds^2 = -dt^2 + t^(2*p1)*dx^2 + t^(2*p2)*dy^2 + t^(2*p3)*dz^2
//
// The parameters p1, p2, p3 can be represented by the Khalatnikov-Lifshitz parameter u:
//   p1 = -u / (1 + u + u^2)
//   p2 = (1 + u) / (1 + u + u^2)
//   p3 = u * (1 + u) / (1 + u + u^2)
//
// The parameters satisfy: p1 + p2 + p3 = 1 and p1^2 + p2^2 + p3^2 = 1
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class kasner : public metric<coordinate_system_type::cartesian, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  kasner()
  {
    calculate_parameters();
  }

  __device__ termination_reason check_termination(const vector_type& position, [[maybe_unused]] const vector_type& direction) const override
  {
    // Check for singularity at t = 0
    if (position[0] <= static_cast<scalar_type>(0))
      return termination_reason::spacetime_breakdown;
    return termination_reason::none;
  }

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto t  = position[0];
    const auto t1 = static_cast<scalar_type>(1) / t;
    const auto t2 = p1 * t1;
    const auto t3 = p2 * t1;
    const auto t4 = p3 * t1;
    const auto t5 = std::pow(t, p1);
    const auto t6 = t5 * t5;
    const auto t9 = std::pow(t, p2);
    const auto t10 = t9 * t9;
    const auto t13 = std::pow(t, p3);
    const auto t14 = t13 * t13;

    christoffel_symbols_type symbols;
    symbols.setZero();
    
    // Non-zero Christoffel symbols
    symbols(0, 1, 1) = t2;
    symbols(0, 2, 2) = t3;
    symbols(0, 3, 3) = t4;
    
    symbols(1, 0, 1) = t2;
    symbols(1, 1, 0) = t6 * p1 * t1;
    
    symbols(2, 0, 2) = t3;
    symbols(2, 2, 0) = t10 * p2 * t1;
    
    symbols(3, 0, 3) = t4;
    symbols(3, 3, 0) = t14 * p3 * t1;
    
    return symbols;
  }

  void set_parameter_u(scalar_type u_value)
  {
    u = u_value;
    calculate_parameters();
  }

  scalar_type get_p1() const { return p1; }
  scalar_type get_p2() const { return p2; }
  scalar_type get_p3() const { return p3; }

private:
  void calculate_parameters()
  {
    const auto denominator = static_cast<scalar_type>(1) + u + u * u;
    p1 = -u / denominator;
    p2 = (static_cast<scalar_type>(1) + u) / denominator;
    p3 = u * (static_cast<scalar_type>(1) + u) / denominator;
  }

  scalar_type u  = static_cast<scalar_type>(0);
  scalar_type p1 = static_cast<scalar_type>(0);
  scalar_type p2 = static_cast<scalar_type>(1);
  scalar_type p3 = static_cast<scalar_type>(0);
};
}
