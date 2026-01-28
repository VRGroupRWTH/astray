#pragma once

#include <cmath>
#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>

namespace ast::metrics
{
// Implementation based on Motion4D metric by Heiko Munz (tauzero7)
// Reference: https://github.com/tauzero7/Motion4D
//
// Plane gravitational wave metric in Cartesian coordinates (t,u,x,y).
// This represents a plane-fronted gravitational wave propagating in the u-direction.
//
// The metric is: ds² = -c² dt² + du² + p²(t,u) dx² + q²(t,u) dy²
// where p and q are wave profile functions that depend on the specific wave form.
//
// NOTE: This implementation requires specification of the wave profile functions p(t,u)
// and q(t,u) and their derivatives. The current implementation uses a simplified form.
//
// Reference: Exact solutions (Stephani et al.)

template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class plane_grav_wave : public metric<coordinate_system_type::cartesian, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  // Wave profile functions - these would be specified based on the particular solution
  __device__ __host__ scalar_type get_p(scalar_type t, scalar_type u) const
  {
    // Simplified profile: p = 1 + amplitude * sin(k * (u - c*t))
    const auto phase = wave_number * (u - consts::speed_of_light * t);
    return static_cast<scalar_type>(1) + amplitude * std::sin(phase);
  }

  __device__ __host__ scalar_type get_q(scalar_type t, scalar_type u) const
  {
    // Simplified profile: q = 1 - amplitude * sin(k * (u - c*t))
    const auto phase = wave_number * (u - consts::speed_of_light * t);
    return static_cast<scalar_type>(1) - amplitude * std::sin(phase);
  }

  __device__ __host__ scalar_type get_dp_du(scalar_type t, scalar_type u) const
  {
    const auto phase = wave_number * (u - consts::speed_of_light * t);
    return amplitude * wave_number * std::cos(phase);
  }

  __device__ __host__ scalar_type get_dq_du(scalar_type t, scalar_type u) const
  {
    const auto phase = wave_number * (u - consts::speed_of_light * t);
    return -amplitude * wave_number * std::cos(phase);
  }

  __device__ termination_reason check_termination(const vector_type& position, const vector_type& direction) const override
  {
    return termination_reason::none;
  }

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto t = position[0];
    const auto u = position[1];
    
    const auto p = get_p(t, u);
    const auto q = get_q(t, u);
    const auto dp = get_dp_du(t, u);
    const auto dq = get_dq_du(t, u);
    
    const auto c = consts::speed_of_light;
    const auto c2 = c * c;
    
    christoffel_symbols_type symbols;
    symbols.setZero();
    
    // Christoffel symbols for plane wave metric
    // These are derived from the metric ds² = -c² dt² + du² + p² dx² + q² dy²
    
    const auto inv_p = static_cast<scalar_type>(1) / p;
    const auto inv_q = static_cast<scalar_type>(1) / q;
    
    // Γ^t_xx = (c/p) * dp/du
    symbols(0, 2, 2) = (c * dp) * inv_p;
    
    // Γ^t_yy = (c/q) * dq/du
    symbols(0, 3, 3) = (c * dq) * inv_q;
    
    // Γ^u_xx = -(1/p) * dp/du
    symbols(1, 2, 2) = -dp * inv_p;
    
    // Γ^u_yy = -(1/q) * dq/du
    symbols(1, 3, 3) = -dq * inv_q;
    
    // Γ^x_tx = Γ^x_xt = (c/p) * dp/du
    symbols(2, 0, 2) = (c * dp) * inv_p;
    symbols(2, 2, 0) = (c * dp) * inv_p;
    
    // Γ^x_ux = Γ^x_xu = (1/p) * dp/du
    symbols(2, 1, 2) = dp * inv_p;
    symbols(2, 2, 1) = dp * inv_p;
    
    // Γ^y_ty = Γ^y_yt = (c/q) * dq/du
    symbols(3, 0, 3) = (c * dq) * inv_q;
    symbols(3, 3, 0) = (c * dq) * inv_q;
    
    // Γ^y_uy = Γ^y_yu = (1/q) * dq/du
    symbols(3, 1, 3) = dq * inv_q;
    symbols(3, 3, 1) = dq * inv_q;
    
    return symbols;
  }
  
  scalar_type amplitude = static_cast<scalar_type>(0.1);     // Wave amplitude
  scalar_type wave_number = static_cast<scalar_type>(1);     // Wave number k
};
}
