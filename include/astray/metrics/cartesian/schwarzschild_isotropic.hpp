#pragma once

#include <cmath>

#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>

namespace ast::metrics
{
// Implementation based on Motion4D metric by Thomas Mueller (tauzero7)
// Reference: https://github.com/tauzero7/Motion4D
//
// Schwarzschild metric in Cartesian isotropic coordinates (t,x,y,z).
//
// The line element is given by:
//   ds^2 = -((1 - ρs/ρ)/(1 + ρs/ρ))^2 c^2 dt^2 + (1 + ρs/ρ)^4 (dx^2 + dy^2 + dz^2)
//
// where ρ^2 = x^2 + y^2 + z^2 and ρs = GM/(2c^2).
//
// This form of the Schwarzschild metric uses isotropic spatial coordinates 
// where the spatial part of the metric is conformally flat.
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class schwarzschild_isotropic : public metric<coordinate_system_type::cartesian, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  __device__ termination_reason check_termination(const vector_type& position, [[maybe_unused]] const vector_type& direction) const override
  {
    const auto rs   = consts::schwarzschild_radius(mass);
    const auto rho_s = static_cast<scalar_type>(0.5) * rs / (consts::speed_of_light_squared);
    const auto rho  = calc_rho(position);
    
    // Check for singularity at rho = rho_s
    if (rho <= rho_s * (static_cast<scalar_type>(1) + consts::epsilon))
      return termination_reason::spacetime_breakdown;
    return termination_reason::none;
  }

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto c  = consts::speed_of_light;
    const auto rs = consts::schwarzschild_radius(mass);
    const auto rho_s = static_cast<scalar_type>(0.5) * rs / (consts::speed_of_light_squared);
    const auto rho = calc_rho(position);
    
    scalar_type drdx, drdy, drdz;
    calc_drho(position, rho, drdx, drdy, drdz);

    const auto t1  = rho * rho;
    const auto t2  = t1 * t1;
    const auto t7  = rho_s * rho_s;
    const auto t13 = t7 * t7;
    const auto t17 = c * c;
    const auto t18 = t2 / (t2 + static_cast<scalar_type>(4) * rho_s * t1 * rho + static_cast<scalar_type>(6) * t7 * t1 + static_cast<scalar_type>(4) * t7 * rho_s * rho + t13) * t17;
    const auto t19 = drdx;
    const auto t20 = rho_s * t19;
    const auto t21 = rho - rho_s;
    const auto t22 = rho + rho_s;
    const auto t23 = t22 * t22;
    const auto t26 = t21 / t23 / t22;
    const auto t30 = drdy;
    const auto t31 = rho_s * t30;
    const auto t35 = drdz;
    const auto t36 = rho_s * t35;
    const auto t40 = static_cast<scalar_type>(1) / t22;
    const auto t42 = t40 / t21;
    const auto t44 = static_cast<scalar_type>(2) * t20 * t42;
    const auto t46 = static_cast<scalar_type>(2) * t31 * t42;
    const auto t48 = static_cast<scalar_type>(2) * t36 * t42;
    const auto t50 = static_cast<scalar_type>(1) / rho * t40;
    const auto t52 = static_cast<scalar_type>(2) * t20 * t50;
    const auto t54 = static_cast<scalar_type>(2) * t31 * t50;
    const auto t56 = static_cast<scalar_type>(2) * t36 * t50;

    christoffel_symbols_type symbols;
    symbols.setZero();
    
    // Non-zero Christoffel symbols
    symbols(0, 0, 1) = static_cast<scalar_type>(2) * t18 * t20 * t26;
    symbols(0, 0, 2) = static_cast<scalar_type>(2) * t18 * t31 * t26;
    symbols(0, 0, 3) = static_cast<scalar_type>(2) * t18 * t36 * t26;
    
    symbols(0, 1, 0) = t44;
    symbols(0, 2, 0) = t46;
    symbols(0, 3, 0) = t48;
    
    symbols(1, 0, 0) = t44;
    symbols(1, 1, 1) = -t52;
    symbols(1, 1, 2) = t54;
    symbols(1, 1, 3) = t56;
    symbols(1, 2, 1) = -t54;
    symbols(1, 2, 2) = -t52;
    symbols(1, 3, 1) = -t56;
    symbols(1, 3, 3) = -t52;
    
    symbols(2, 0, 0) = t46;
    symbols(2, 1, 1) = -t54;
    symbols(2, 1, 2) = -t52;
    symbols(2, 2, 1) = t52;
    symbols(2, 2, 2) = -t54;
    symbols(2, 2, 3) = t56;
    symbols(2, 3, 2) = -t56;
    symbols(2, 3, 3) = -t54;
    
    symbols(3, 0, 0) = t48;
    symbols(3, 1, 1) = -t56;
    symbols(3, 1, 3) = -t52;
    symbols(3, 2, 2) = -t56;
    symbols(3, 2, 3) = -t54;
    symbols(3, 3, 1) = t52;
    symbols(3, 3, 2) = t54;
    symbols(3, 3, 3) = -t56;
    
    return symbols;
  }

private:
  __device__ scalar_type calc_rho(const vector_type& position) const
  {
    const auto x = position[1];
    const auto y = position[2];
    const auto z = position[3];
    return std::sqrt(x * x + y * y + z * z);
  }

  __device__ void calc_drho(const vector_type& position, scalar_type rho, scalar_type& drdx, scalar_type& drdy, scalar_type& drdz) const
  {
    const auto x = position[1];
    const auto y = position[2];
    const auto z = position[3];
    const auto inv_rho = static_cast<scalar_type>(1) / rho;
    drdx = x * inv_rho;
    drdy = y * inv_rho;
    drdz = z * inv_rho;
  }

public:
  scalar_type mass = static_cast<scalar_type>(1);
};
}
