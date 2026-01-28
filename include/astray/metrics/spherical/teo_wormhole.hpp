#pragma once

#include <cmath>

#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>

namespace ast::metrics
{
// Implementation based on Motion4D metric by Thomas Mueller (tauzero7)
// Reference: https://github.com/tauzero7/Motion4D
//
// Teo metric of an axisymmetric rotating wormhole in spherical coordinates (t,l,theta,phi).
//
// The line element is given by:
//   ds² = -N(l,θ) c² dt² + dl² + r(l)² K(l,θ)² [dθ² + sin²(θ)(dφ - ω(l,θ) c dt)²]
//
// The standard potentials are:
//   N(l,θ) = 1, K(l,θ) = 1, r(l) = √(l² + b₀²), ω(l,θ) = b₀²/(2(l² + b₀²)^(3/2))
//
// References:
//   Edward Teo, "Rotating traversable wormholes", Phys. Rev. D 58, 024014 (1998).
//   Oliver Fechtig, "Physikalische Aspekte und Visualisierung von stationaeren Wurmloechern"
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class teo_wormhole : public metric<coordinate_system_type::spherical, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto l     = position[1];
    const auto theta = position[2];
    const auto c     = consts::speed_of_light;
    
    // Calculate potentials and their derivatives
    const auto N  = static_cast<scalar_type>(1);
    const auto dN = static_cast<scalar_type>(0);
    const auto K  = static_cast<scalar_type>(1);
    const auto dK = static_cast<scalar_type>(0);
    
    const auto b0_sq = b0 * b0;
    const auto l_sq  = l * l;
    const auto r_sq_sum = b0_sq + l_sq;
    
    const auto r  = std::sqrt(r_sq_sum);
    const auto dr = l / r;
    const auto d2r = b0_sq / std::pow(r_sq_sum, static_cast<scalar_type>(1.5));
    
    const auto omega  = static_cast<scalar_type>(0.5) * b0_sq * std::pow(r, static_cast<scalar_type>(-3));
    const auto domega = static_cast<scalar_type>(-1.5) * b0_sq * l * std::pow(r_sq_sum, static_cast<scalar_type>(-2.5));
    
    // Precompute common terms
    const auto st = std::sin(theta);
    const auto ct = std::cos(theta);
    const auto st_sq = st * st;
    const auto c_sq = c * c;
    
    const auto t1 = N;
    const auto t2 = c_sq;
    const auto t4 = dN;
    const auto t6 = r;
    const auto t7 = K;
    const auto t8 = t7 * t7;
    const auto t9 = t6 * t8;
    const auto t10 = st;
    const auto t11 = t10 * t10;
    const auto t12 = t9 * t11;
    const auto t13 = omega;
    const auto t14 = t13 * t13;
    const auto t15 = t14 * t2;
    const auto t16 = dr;
    const auto t19 = t6 * t6;
    const auto t20 = t19 * t7;
    const auto t21 = t20 * t11;
    const auto t22 = dK;
    const auto t25 = t19 * t8;
    const auto t26 = t25 * t11;
    const auto t28 = domega;
    const auto t33 = ct;
    const auto t42 = t1 * t1;
    const auto t43 = static_cast<scalar_type>(1) / t42;
    const auto t45 = (static_cast<scalar_type>(2) * t1 * t4 - t25 * t11 * t13 * t28) * t43 / static_cast<scalar_type>(2);
    const auto t46 = static_cast<scalar_type>(1) / t6;
    const auto t48 = static_cast<scalar_type>(1) / t7;
    const auto t54 = t19 * t6;
    const auto t58 = t8 * t7 * t11 * t28;
    const auto t60 = t42 * t7;
    const auto t64 = t42 * t6;
    const auto t74 = c * t46 * t48
        * (static_cast<scalar_type>(2) * t13 * t6 * t7 * t1 * t4 - t14 * t54 * t58 - static_cast<scalar_type>(2) * t60 * t13 * t16 
           - static_cast<scalar_type>(2) * t64 * t13 * t22 - t64 * t7 * t28)
        * t43 / static_cast<scalar_type>(2);
    const auto t75 = t13 * c;
    const auto t77 = t33 / t10;
    const auto t78 = t75 * t77;
    const auto t87 = t12 * t75 * t16 + t21 * t75 * t22 + t25 * t11 * t28 * c / static_cast<scalar_type>(2);
    const auto t90 = t10 * t13 * c * t33;
    const auto t91 = t46 * t48;
    const auto t95 = t91 * (t7 * t16 + t6 * t22);
    const auto t100 = t26 * t28 / c * t43 / static_cast<scalar_type>(2);
    const auto t110 = t91 * (t13 * t54 * t58 + static_cast<scalar_type>(2) * t60 * t16 + static_cast<scalar_type>(2) * t64 * t22) * t43 / static_cast<scalar_type>(2);

    christoffel_symbols_type symbols;
    symbols.setZero();
    
    // Non-zero Christoffel symbols
    symbols(0, 0, 1) = t1 * t2 * t4 - t12 * t15 * t16 - t21 * t15 * t22 - t26 * t13 * t2 * t28;
    symbols(0, 0, 2) = -t10 * t14 * t2 * t33;
    
    symbols(0, 1, 0) = t45;
    symbols(0, 1, 3) = t74;
    
    symbols(0, 2, 3) = -t78;
    
    symbols(0, 3, 1) = t87;
    symbols(0, 3, 2) = t90;
    
    symbols(1, 0, 0) = t45;
    symbols(1, 0, 3) = t74;
    
    symbols(1, 2, 2) = t95;
    
    symbols(1, 3, 0) = t100;
    symbols(1, 3, 3) = t110;
    
    symbols(2, 0, 3) = -t78;
    
    symbols(2, 1, 2) = t95;
    
    symbols(2, 2, 1) = -t9 * t16 - t20 * t22;
    
    symbols(2, 3, 3) = t77;
    
    symbols(3, 0, 1) = t87;
    symbols(3, 0, 2) = t90;
    
    symbols(3, 1, 0) = t100;
    symbols(3, 1, 3) = t110;
    
    symbols(3, 2, 3) = t77;
    
    symbols(3, 3, 1) = -t9 * t11 * t16 - t20 * t11 * t22;
    symbols(3, 3, 2) = -t10 * t33;
    
    return symbols;
  }

  // Wormhole throat radius parameter
  scalar_type b0 = static_cast<scalar_type>(1);
};
}
