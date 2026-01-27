#pragma once

#include <cmath>

#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>

namespace ast::metrics
{
// Implementation based on Motion4D metric by Thomas Mueller (tauzero7)
// Reference: https://github.com/tauzero7/Motion4D
//
// Rotating extreme Reissner-Nordstrom dihole metric in Cartesian coordinates (t,x,y,z).
//
// The metric is given by:
//   g_tt = -1/U², g_xx = g_yy = g_zz = U²
// where U = 1 + M1/r1 + M2/r2
//
// The two holes are located at positions that rotate with angular velocity ω:
//   r1 = sqrt(x² + (y + sin(ωt))² + (z - cos(ωt))²)
//   r2 = sqrt(x² + (y - sin(ωt))² + (z + cos(ωt))²)
//
// Parameters:
//   mass1, mass2 - masses of the two holes
//   omega - angular velocity of rotation
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class rot_dihole : public metric<coordinate_system_type::cartesian, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  __device__ termination_reason check_termination(const vector_type& position, [[maybe_unused]] const vector_type& direction) const override
  {
    const auto t = position[0];
    const auto x = position[1];
    const auto y = position[2];
    const auto z = position[3];
    
    const auto omega = angular_velocity;
    const auto sin_omega_t = std::sin(omega * t);
    const auto cos_omega_t = std::cos(omega * t);
    const auto r1 = std::sqrt(x * x + (y + sin_omega_t) * (y + sin_omega_t) + (z - cos_omega_t) * (z - cos_omega_t));
    const auto r2 = std::sqrt(x * x + (y - sin_omega_t) * (y - sin_omega_t) + (z + cos_omega_t) * (z + cos_omega_t));
    
    // Check for singularities at the hole locations
    if (r1 < consts::epsilon || r2 < consts::epsilon)
      return termination_reason::spacetime_breakdown;
    
    return termination_reason::none;
  }

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto t = position[0];
    const auto x = position[1];
    const auto y = position[2];
    const auto z = position[3];
    
    const auto M1 = mass1;
    const auto M2 = mass2;
    const auto omega = angular_velocity;
    
    // Calculate r1 and r2
    const auto sin_omega_t = std::sin(omega * t);
    const auto cos_omega_t = std::cos(omega * t);
    const auto r1 = std::sqrt(x * x + (y + sin_omega_t) * (y + sin_omega_t) + (z - cos_omega_t) * (z - cos_omega_t));
    const auto r2 = std::sqrt(x * x + (y - sin_omega_t) * (y - sin_omega_t) + (z + cos_omega_t) * (z + cos_omega_t));
    
    // Calculate U
    const auto U = static_cast<scalar_type>(1) + M1 / r1 + M2 / r2;
    
    // Calculate derivatives of U
    const auto r1_3 = r1 * r1 * r1;
    const auto r2_3 = r2 * r2 * r2;
    const auto dUdt = -M1 * omega * ((y + sin_omega_t) * cos_omega_t + (z - cos_omega_t) * sin_omega_t) / r1_3
                    + M2 * omega * ((y - sin_omega_t) * cos_omega_t + (z + cos_omega_t) * sin_omega_t) / r2_3;
    const auto dUdx = -M1 * x / r1_3 - M2 * x / r2_3;
    const auto dUdy = -M1 * (y + sin_omega_t) / r1_3 - M2 * (y - sin_omega_t) / r2_3;
    const auto dUdz = -M1 * (z - cos_omega_t) / r1_3 - M2 * (z + cos_omega_t) / r2_3;
    
    const auto U_inv = static_cast<scalar_type>(1) / U;
    const auto U_3 = U * U * U;
    const auto U_5 = U_3 * U * U;
    
    christoffel_symbols_type symbols;
    symbols.setZero();
    
    // Non-zero Christoffel symbols
    symbols(0, 0, 0) = -dUdt * U_inv;
    symbols(0, 0, 1) = -dUdx / U_5;
    symbols(0, 0, 2) = -dUdy / U_5;
    symbols(0, 0, 3) = -dUdz / U_5;
    symbols(0, 1, 0) = -dUdx * U_inv;
    symbols(0, 1, 1) = dUdt * U_inv;
    symbols(0, 2, 0) = -dUdy * U_inv;
    symbols(0, 2, 2) = dUdt * U_inv;
    symbols(0, 3, 0) = -dUdz * U_inv;
    symbols(0, 3, 3) = dUdt * U_inv;
    
    symbols(1, 0, 0) = -dUdx * U_inv;
    symbols(1, 0, 1) = dUdt * U_inv;
    symbols(1, 1, 0) = U_3 * dUdt;
    symbols(1, 1, 1) = dUdx * U_inv;
    symbols(1, 1, 2) = -dUdy * U_inv;
    symbols(1, 1, 3) = -dUdz * U_inv;
    symbols(1, 2, 1) = dUdy * U_inv;
    symbols(1, 2, 2) = dUdx * U_inv;
    symbols(1, 3, 1) = dUdz * U_inv;
    symbols(1, 3, 3) = dUdx * U_inv;
    
    symbols(2, 0, 0) = -dUdy * U_inv;
    symbols(2, 0, 2) = dUdt * U_inv;
    symbols(2, 1, 1) = dUdy * U_inv;
    symbols(2, 1, 2) = dUdx * U_inv;
    symbols(2, 2, 0) = U_3 * dUdt;
    symbols(2, 2, 1) = -dUdx * U_inv;
    symbols(2, 2, 2) = dUdy * U_inv;
    symbols(2, 2, 3) = -dUdz * U_inv;
    symbols(2, 3, 2) = dUdz * U_inv;
    symbols(2, 3, 3) = dUdy * U_inv;
    
    symbols(3, 0, 0) = -dUdz * U_inv;
    symbols(3, 0, 3) = dUdt * U_inv;
    symbols(3, 1, 1) = dUdz * U_inv;
    symbols(3, 1, 3) = dUdx * U_inv;
    symbols(3, 2, 2) = dUdz * U_inv;
    symbols(3, 2, 3) = dUdy * U_inv;
    symbols(3, 3, 0) = U_3 * dUdt;
    symbols(3, 3, 1) = -dUdx * U_inv;
    symbols(3, 3, 2) = -dUdy * U_inv;
    symbols(3, 3, 3) = dUdz * U_inv;
    
    return symbols;
  }

  scalar_type mass1 = static_cast<scalar_type>(1);
  scalar_type mass2 = static_cast<scalar_type>(1);
  scalar_type angular_velocity = static_cast<scalar_type>(0);
};
}
