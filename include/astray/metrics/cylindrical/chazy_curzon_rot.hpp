#pragma once

#include <cmath>

#include <astray/core/metric.hpp>

namespace ast::metrics
{
// Implementation based on Motion4D metric by Thomas Mueller (tauzero7)
// Reference: https://github.com/tauzero7/Motion4D
//
// Chazy-Curzon rotating metric in Weyl cylindrical coordinates (t,ρ,φ,z).
//
// The line element is given by:
//   ds² = e^(-2U)[e^(2k)(dρ² + dz²) + ρ²dφ²] - e^(2U)(dt + A dφ)²
//
// where:
//   e^(-2U) = cosh(2m/r) - p sinh(2m/r)
//   2k = -m²ρ²/r⁴
//   A = 2qmz/r
//   r² = ρ² + z²
//
// The parameters p, q are related via p² + q² = 1.
// This represents a rotating generalization of the Curzon metric.
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class chazy_curzon_rot : public metric<coordinate_system_type::cylindrical, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto rho = position[1];
    
    scalar_type em2U, k, A;
    calc_UkA(position, em2U, k, A);
    const auto U = static_cast<scalar_type>(-0.5) * std::log(em2U);
    
    scalar_type dUdrho, dUdz, dkdrho, dkdz, dAdrho, dAdz;
    calc_DUka(position, dUdrho, dUdz, dkdrho, dkdz, dAdrho, dAdz);

    const auto t1  = k;
    const auto t2  = std::exp(t1);
    const auto t3  = t2 * t2;
    const auto t4  = static_cast<scalar_type>(1) / t3;
    const auto t5  = U;
    const auto t6  = std::exp(t5);
    const auto t7  = t6 * t6;
    const auto t8  = t7 * t7;
    const auto t9  = t4 * t8;
    const auto t10 = dUdrho;
    const auto t12 = dUdz;
    const auto t14 = rho * rho;
    const auto t15 = t14 * t10;
    const auto t16 = static_cast<scalar_type>(2) * t15;
    const auto t17 = A;
    const auto t18 = t8 * t17;
    const auto t19 = dAdrho;
    const auto t20 = t18 * t19;
    const auto t22 = static_cast<scalar_type>(1) / t14;
    const auto t24 = (t16 + t20) * t22 / static_cast<scalar_type>(2);
    const auto t27 = t8 * t19 * t22 / static_cast<scalar_type>(2);
    const auto t33 = t8 * (static_cast<scalar_type>(2) * t17 * t10 + t19) * t4 / static_cast<scalar_type>(2);
    const auto t36 = dAdz;
    const auto t40 = t8 * (static_cast<scalar_type>(2) * t17 * t12 + t36) * t4 / static_cast<scalar_type>(2);
    const auto t41 = t14 * t12;
    const auto t43 = t18 * t36;
    const auto t46 = (static_cast<scalar_type>(2) * t41 + t43) * t22 / static_cast<scalar_type>(2);
    const auto t49 = t8 * t36 * t22 / static_cast<scalar_type>(2);
    const auto t50 = dkdrho;
    const auto t51 = -t10 + t50;
    const auto t52 = dkdz;
    const auto t53 = t12 - t52;
    const auto t54 = t14 * t17;
    const auto t58 = t17 * t17;
    const auto t59 = t8 * t58;
    const auto t65 = (static_cast<scalar_type>(4) * t54 * t10 + t14 * t19 + t59 * t19 - static_cast<scalar_type>(2) * t17 * rho) * t22 / static_cast<scalar_type>(2);
    const auto t69 = (t20 + t16 - static_cast<scalar_type>(2) * rho) * t22 / static_cast<scalar_type>(2);
    const auto t82 = (static_cast<scalar_type>(4) * t54 * t12 + t14 * t36 + t59 * t36) * t22 / static_cast<scalar_type>(2);

    christoffel_symbols_type symbols;
    symbols.setZero();
    
    // Non-zero Christoffel symbols
    symbols(0, 0, 1) = t9 * t10;
    symbols(0, 0, 3) = t9 * t12;
    symbols(0, 1, 0) = t24;
    symbols(0, 1, 2) = -t27;
    symbols(0, 2, 1) = t33;
    symbols(0, 2, 3) = t40;
    symbols(0, 3, 0) = t46;
    symbols(0, 3, 2) = -t49;
    
    symbols(1, 0, 0) = t24;
    symbols(1, 0, 2) = -t27;
    symbols(1, 1, 1) = t51;
    symbols(1, 1, 3) = t53;
    symbols(1, 2, 0) = t65;
    symbols(1, 2, 2) = -t69;
    symbols(1, 3, 1) = -t53;
    symbols(1, 3, 3) = t51;
    
    symbols(2, 0, 1) = t33;
    symbols(2, 0, 3) = t40;
    symbols(2, 1, 0) = t65;
    symbols(2, 1, 2) = -t69;
    symbols(2, 2, 1) = t4 * (t15 - rho + t59 * t10 + t20);
    symbols(2, 2, 3) = t4 * (t41 + t59 * t12 + t43);
    symbols(2, 3, 0) = t82;
    symbols(2, 3, 2) = -t46;
    
    symbols(3, 0, 0) = t46;
    symbols(3, 0, 2) = -t49;
    symbols(3, 1, 1) = -t53;
    symbols(3, 1, 3) = t51;
    symbols(3, 2, 0) = t82;
    symbols(3, 2, 2) = -t46;
    symbols(3, 3, 1) = -t51;
    symbols(3, 3, 3) = -t53;
    
    return symbols;
  }

  void set_parameters(scalar_type mass_val, scalar_type p_val)
  {
    mass = mass_val;
    p = p_val;
    q = std::sqrt(static_cast<scalar_type>(1) - p * p);
  }

private:
  __device__ void calc_UkA(const vector_type& position, scalar_type& em2U, scalar_type& k, scalar_type& A) const
  {
    const auto rho = position[1];
    const auto z = position[3];
    const auto r2 = rho * rho + z * z;
    const auto r = std::sqrt(r2);

    em2U = std::cosh(static_cast<scalar_type>(2) * mass / r) - p * std::sinh(static_cast<scalar_type>(2) * mass / r);
    k = static_cast<scalar_type>(-0.5) * mass * mass * rho * rho / (r2 * r2);
    A = static_cast<scalar_type>(2) * q * mass * z / r;
  }

  __device__ void calc_DUka(const vector_type& position, scalar_type& dUdrho, scalar_type& dUdz, 
                             scalar_type& dkdrho, scalar_type& dkdz, scalar_type& dAdrho, scalar_type& dAdz) const
  {
    const auto rho = position[1];
    const auto z = position[3];
    const auto r2 = rho * rho + z * z;
    const auto r = std::sqrt(r2);

    const auto ch = std::cosh(static_cast<scalar_type>(2) * mass / r);
    const auto sh = std::sinh(static_cast<scalar_type>(2) * mass / r);
    const auto w1 = r * r * r;

    dUdrho = mass * rho * (sh - p * ch) / (w1 * (ch - p * sh));
    dUdz = mass * z * (sh - p * ch) / (w1 * (ch - p * sh));

    dkdrho = mass * mass * rho * (rho * rho - z * z) / (r2 * r2 * r2);
    dkdz = static_cast<scalar_type>(2) * mass * mass * rho * rho * z / (r2 * r2 * r2);

    dAdrho = static_cast<scalar_type>(-2) * mass * q * rho * z / (r * r * r);
    dAdz = static_cast<scalar_type>(2) * mass * q * rho * rho / (r * r * r);
  }

  scalar_type mass = static_cast<scalar_type>(1);
  scalar_type p = static_cast<scalar_type>(0.5);
  scalar_type q = static_cast<scalar_type>(0.866025403784439); // sqrt(1 - 0.5^2)
};
}
