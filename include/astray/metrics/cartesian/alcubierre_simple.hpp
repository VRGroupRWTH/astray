#pragma once

#include <cmath>

#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>

namespace ast::metrics
{
// Implementation based on Motion4D metric by Thomas Mueller (tauzero7)
// Reference: https://github.com/tauzero7/Motion4D
//
// Alcubierre warp drive metric (simplified version) in Cartesian coordinates (t,x,y,z).
//
// The line element is given by:
//   ds² = -c²dt² + [dx - vs f(rs) dt]² + dy² + dz²
//
// where:
//   vs = velocity of warp bubble
//   rs = √((x - vs*t)² + y² + z²) = distance from bubble center
//   f(rs) = shape function that creates the warp bubble
//
// For the simplified version with smooth transition:
//   f(rs) = { 1,                    rs ≤ R1
//           { (1 - w²)²,            R1 < rs < R2  where w = (rs - R1)/(R2 - R1)
//           { 0,                    rs ≥ R2
//
// where R1 = R - dR/2, R2 = R + dR/2
//
// Reference: M. Alcubierre, Classical Quantum Gravity 11, L73-L77 (1994)
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class alcubierre_simple : public metric<coordinate_system_type::cartesian, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto c = consts::speed_of_light;
    const auto f = calc_f(position);
    
    scalar_type ft, fx, fy, fz;
    calc_df(position, ft, fx, fy, fz);

    const auto vs2 = vs * vs;
    const auto vs3 = vs2 * vs;
    const auto c2 = c * c;
    const auto f2 = f * f;
    const auto f3 = f2 * f;

    christoffel_symbols_type symbols;
    symbols.setZero();
    
    // Non-zero Christoffel symbols
    symbols(0, 0, 0) = vs3 * f2 * fx / c2;
    symbols(0, 0, 1) = vs * (-c2 * vs * f * fx - c2 * ft + vs3 * f3 * fx) / c2;
    symbols(0, 0, 2) = -vs2 * f * fy;
    symbols(0, 0, 3) = -vs2 * f * fz;
    
    symbols(0, 1, 0) = -vs2 * f * fx / c2;
    symbols(0, 1, 1) = -vs3 * f2 * fx / c2;
    symbols(0, 1, 2) = vs * fy / static_cast<scalar_type>(2);
    symbols(0, 1, 3) = vs * fz / static_cast<scalar_type>(2);
    
    symbols(0, 2, 0) = -vs2 * f * fy / (static_cast<scalar_type>(2) * c2);
    symbols(0, 2, 1) = vs * (-c2 - vs2 * f2) * fy / (static_cast<scalar_type>(2) * c2);
    
    symbols(0, 3, 0) = -vs2 * f * fz / (static_cast<scalar_type>(2) * c2);
    symbols(0, 3, 1) = vs * (-c2 - vs2 * f2) * fz / (static_cast<scalar_type>(2) * c2);
    
    symbols(1, 0, 0) = -vs2 * f * fx / c2;
    symbols(1, 0, 1) = -vs3 * f2 * fx / c2;
    symbols(1, 0, 2) = vs * fy / static_cast<scalar_type>(2);
    symbols(1, 0, 3) = vs * fz / static_cast<scalar_type>(2);
    
    symbols(1, 1, 0) = vs * fx / c2;
    symbols(1, 1, 1) = vs2 * f * fx / c2;
    
    symbols(1, 2, 0) = vs * fy / (static_cast<scalar_type>(2) * c2);
    symbols(1, 2, 1) = vs2 * f * fy / (static_cast<scalar_type>(2) * c2);
    
    symbols(1, 3, 0) = vs * fz / (static_cast<scalar_type>(2) * c2);
    symbols(1, 3, 1) = vs2 * f * fz / (static_cast<scalar_type>(2) * c2);
    
    symbols(2, 0, 0) = -vs2 * f * fy / (static_cast<scalar_type>(2) * c2);
    symbols(2, 0, 1) = vs * (-c2 - vs2 * f2) * fy / (static_cast<scalar_type>(2) * c2);
    
    symbols(2, 1, 0) = vs * fy / (static_cast<scalar_type>(2) * c2);
    symbols(2, 1, 1) = vs2 * f * fy / (static_cast<scalar_type>(2) * c2);
    
    symbols(3, 0, 0) = -vs2 * f * fz / (static_cast<scalar_type>(2) * c2);
    symbols(3, 0, 1) = vs * (-c2 - vs2 * f2) * fz / (static_cast<scalar_type>(2) * c2);
    
    symbols(3, 1, 0) = vs * fz / (static_cast<scalar_type>(2) * c2);
    symbols(3, 1, 1) = vs2 * f * fz / (static_cast<scalar_type>(2) * c2);
    
    return symbols;
  }

  void set_parameters(scalar_type R_val, scalar_type dR_val, scalar_type vs_val)
  {
    R = R_val;
    dR = dR_val;
    vs = vs_val;
  }

private:
  __device__ scalar_type calc_rs(const vector_type& position) const
  {
    const auto t = position[0];
    const auto x = position[1] - vs * t;
    const auto y = position[2];
    const auto z = position[3];
    return std::sqrt(x * x + y * y + z * z);
  }

  __device__ scalar_type calc_f(const vector_type& position) const
  {
    const auto rs = calc_rs(position);
    const auto R1 = R - static_cast<scalar_type>(0.5) * dR;
    const auto R2 = R + static_cast<scalar_type>(0.5) * dR;
    
    if (rs <= R1)
      return static_cast<scalar_type>(1);
    
    if (rs >= R1 && rs < R2)
    {
      const auto w = (rs - R1) / (R2 - R1);
      const auto y = static_cast<scalar_type>(1) - w * w;
      return y * y;
    }
    
    return static_cast<scalar_type>(0);
  }

  __device__ void calc_df(const vector_type& position, scalar_type& ft, scalar_type& fx, scalar_type& fy, scalar_type& fz) const
  {
    const auto rs = calc_rs(position);
    const auto t = position[0];
    const auto x = position[1] - vs * t;
    const auto y = position[2];
    const auto z = position[3];
    
    const auto R1 = R - static_cast<scalar_type>(0.5) * dR;
    const auto R2 = R + static_cast<scalar_type>(0.5) * dR;
    
    scalar_type dfdr = static_cast<scalar_type>(0);
    
    if (rs > R1 && rs < R2)
    {
      const auto w = (rs - R1) / (R2 - R1);
      const auto y_val = static_cast<scalar_type>(1) - w * w;
      dfdr = static_cast<scalar_type>(-4) * w * y_val / (R2 - R1);
    }
    
    if (rs > static_cast<scalar_type>(0))
    {
      const auto df = dfdr / rs;
      ft = -vs * x * df;
      fx = x * df;
      fy = y * df;
      fz = z * df;
    }
    else
    {
      ft = fx = fy = fz = static_cast<scalar_type>(0);
    }
  }

  scalar_type R  = static_cast<scalar_type>(1);
  scalar_type dR = static_cast<scalar_type>(0.1);
  scalar_type vs = static_cast<scalar_type>(1);
};
}
