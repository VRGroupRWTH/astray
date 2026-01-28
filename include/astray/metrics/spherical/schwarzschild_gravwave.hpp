#pragma once

#include <cmath>
#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>

namespace ast::metrics
{
// Implementation based on Motion4D metric by Thomas Mueller (tauzero7)
// Reference: https://github.com/tauzero7/Motion4D
//
// Schwarzschild metric with gravitational wave perturbations in spherical coordinates.
//
// Line element (with perturbations):
//   ds^2 = -(1-rs/r + ε*h_tt) dt^2 + (1/(1-rs/r) + ε*h_rr) dr^2 
//          + (r^2 + ε*h_θθ) dθ^2 + (r^2 sin²θ + ε*h_φφ) dφ^2
//
// where ε is the perturbation amplitude, σ is the wave frequency,
// and h_αβ are the perturbation components that depend on Legendre polynomials.
//
// Note: This implementation supports l=0 and l=1 Legendre polynomials only.
// For higher-order modes, the Legendre polynomial calculations need to be extended.

template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class schwarzschild_gravwave : public metric<coordinate_system_type::spherical, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

private:
  // Helper function to calculate Legendre polynomials and their derivatives
  // Only supports l=0 and l=1
  __device__ void calc_legendre(int l, scalar_type theta, 
                                scalar_type& Pl, scalar_type& Plt, 
                                scalar_type& Pltt, scalar_type& Plttt) const
  {
    if (l == 0) {
      Pl = static_cast<scalar_type>(1);
      Plt = static_cast<scalar_type>(0);
      Pltt = static_cast<scalar_type>(0);
      Plttt = static_cast<scalar_type>(0);
    }
    else if (l == 1) {
      const auto ct = std::cos(theta);
      const auto st = std::sin(theta);
      Pl = ct;                    // P_1(cos θ) = cos θ
      Plt = -st;                  // dP_1/dθ = -sin θ
      Pltt = -ct;                 // d²P_1/dθ² = -cos θ
      Plttt = st;                 // d³P_1/dθ³ = sin θ
    }
    else {
      // Higher-order Legendre polynomials not implemented
      Pl = static_cast<scalar_type>(0);
      Plt = static_cast<scalar_type>(0);
      Pltt = static_cast<scalar_type>(0);
      Plttt = static_cast<scalar_type>(0);
    }
  }

  // Calculate perturbations h_αβ
  __device__ void calc_perturbations(const vector_type& position,
                                     scalar_type& htt, scalar_type& hrr, 
                                     scalar_type& hee, scalar_type& hpp) const
  {
    const auto t = position[0];
    const auto r = position[1];
    const auto theta = position[2];
    
    const auto rs = consts::schwarzschild_radius(mass);
    const auto M = static_cast<scalar_type>(0.5) * rs;
    const auto r2 = r * r;
    const auto r4 = r2 * r2;
    
    const auto f = static_cast<scalar_type>(1) - rs / r;
    const auto p = M - (M * M + sigma * sigma * r4) / (r - static_cast<scalar_type>(2) * M);
    const auto q = std::sqrt(f) / r2;
    
    const auto X = p * q;
    const auto Y = static_cast<scalar_type>(3) * M * q;
    const auto Z = (r - static_cast<scalar_type>(3) * M) * q;
    const auto W = r * q;
    
    scalar_type Pl, Plt, Pltt, Plttt;
    calc_legendre(l, theta, Pl, Plt, Pltt, Plttt);
    
    const auto cst = std::cos(sigma * t);
    const auto sth = std::sin(theta);
    
    htt = -f * X * Pl * cst;
    hrr = Y / f * Pl * cst;
    hee = r2 * (Z * Pl + W * Pltt) * cst;
    hpp = r2 * sth * sth * (Z * Pl + W * Plt * std::cos(theta) / sth) * cst;
  }

  // Calculate perturbations and their derivatives
  __device__ void calc_perturbations_and_diffs(const vector_type& position,
                                                scalar_type& htt, scalar_type& hrr, scalar_type& hee, scalar_type& hpp,
                                                scalar_type& htt_t, scalar_type& htt_r, scalar_type& htt_theta,
                                                scalar_type& hrr_t, scalar_type& hrr_r, scalar_type& hrr_theta,
                                                scalar_type& hee_t, scalar_type& hee_r, scalar_type& hee_theta,
                                                scalar_type& hpp_t, scalar_type& hpp_r, scalar_type& hpp_theta) const
  {
    const auto t = position[0];
    const auto r = position[1];
    const auto theta = position[2];
    
    const auto rs = consts::schwarzschild_radius(mass);
    const auto M = static_cast<scalar_type>(0.5) * rs;
    const auto r2 = r * r;
    const auto r3 = r2 * r;
    const auto r4 = r2 * r2;
    
    const auto f = static_cast<scalar_type>(1) - rs / r;
    const auto df = rs / r2;
    const auto sf = std::sqrt(f);
    const auto p = M - (M * M + sigma * sigma * r4) / (r - static_cast<scalar_type>(2) * M);
    const auto q = sf / r2;
    
    const auto X = p * q;
    const auto Y = static_cast<scalar_type>(3) * M * q;
    const auto Z = (r - static_cast<scalar_type>(3) * M) * q;
    const auto W = r * q;
    
    scalar_type Pl, Plt, Pltt, Plttt;
    calc_legendre(l, theta, Pl, Plt, Pltt, Plttt);
    
    const auto cst = std::cos(sigma * t);
    const auto sst = std::sin(sigma * t);
    const auto sth = std::sin(theta);
    const auto cth = std::cos(theta);
    const auto cot = cth / sth;
    
    htt = -f * X * Pl * cst;
    hrr = Y / f * Pl * cst;
    hee = r2 * (Z * Pl + W * Pltt) * cst;
    hpp = r2 * sth * sth * (Z * Pl + W * Plt * cot) * cst;
    
    const auto dp = -(sigma * sigma * r3 * (static_cast<scalar_type>(3) * r - static_cast<scalar_type>(4) * rs) - M * M) / ((r - rs) * (r - rs));
    const auto dq = (static_cast<scalar_type>(0.5) * rs / sf - static_cast<scalar_type>(2) * sf * r) / r4;
    const auto DX = dp * q + p * dq;
    const auto DY = static_cast<scalar_type>(3) * M * dq;
    const auto DZ = q + (r - static_cast<scalar_type>(3) * M) * dq;
    const auto DW = q + r * dq;
    
    htt_t = f * X * Pl * sst * sigma;
    htt_r = -(df * X + f * DX) * Pl * cst;
    htt_theta = -f * X * Plt * cst;
    
    hrr_t = -Y * Pl * sst * sigma / f;
    hrr_r = -Y * Pl * cst * df / (f * f) + DY * Pl * cst / f;
    hrr_theta = Y * Plt * cst / f;
    
    hee_t = -r2 * (Z * Pl + W * Pltt) * sst * sigma;
    hee_r = static_cast<scalar_type>(2) * r * (Z * Pl + W * Pltt) * cst + r2 * (DZ * Pl + DW * Pltt) * cst;
    hee_theta = r2 * (Z * Plt + W * Plttt) * cst;
    
    hpp_t = -r2 * sth * sth * (Z * Pl + W * Plt * cot) * sst * sigma;
    hpp_r = static_cast<scalar_type>(2) * r * sth * sth * (Z * Pl + W * Plt * cot) * cst + r2 * sth * sth * (DZ * Pl + DW * Plt * cot) * cst;
    hpp_theta = static_cast<scalar_type>(2) * r2 * sth * (Z * Pl + W * Plt * cot) * cst * cth
              + r2 * sth * sth * (Z * Plt + W * Pltt * cot + W * Plt * (-static_cast<scalar_type>(1) - cot * cot)) * cst;
  }

public:
  __device__ termination_reason check_termination(const vector_type& position, const vector_type& direction) const override
  {
    const auto rs = consts::schwarzschild_radius(mass);
    if (position[1] < static_cast<scalar_type>(0) || 
        static_cast<scalar_type>(std::pow(position[1], 2)) <= (static_cast<scalar_type>(1) + consts::epsilon) * static_cast<scalar_type>(std::pow(rs, 2)))
      return termination_reason::spacetime_breakdown;
    return termination_reason::none;
  }

  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto r = position[1];
    const auto theta = position[2];
    const auto rs = consts::schwarzschild_radius(mass);
    const auto f = static_cast<scalar_type>(1) - rs / r;
    const auto df = rs / (r * r);
    
    // Calculate perturbations and their derivatives
    scalar_type htt, hrr, hee, hpp;
    scalar_type htt_t, htt_r, htt_theta;
    scalar_type hrr_t, hrr_r, hrr_theta;
    scalar_type hee_t, hee_r, hee_theta;
    scalar_type hpp_t, hpp_r, hpp_theta;
    
    calc_perturbations_and_diffs(position, 
                                 htt, hrr, hee, hpp,
                                 htt_t, htt_r, htt_theta,
                                 hrr_t, hrr_r, hrr_theta,
                                 hee_t, hee_r, hee_theta,
                                 hpp_t, hpp_r, hpp_theta);
    
    const auto c = consts::speed_of_light;
    const auto c2 = c * c;
    
    // Calculate inverse metric factors (these appear in Christoffel symbols)
    const auto t1 = f;
    const auto t2 = htt;
    const auto t5 = static_cast<scalar_type>(1) / (t1 - epsilon * t2);
    const auto t6 = t5 * epsilon;
    const auto t7 = htt_t;
    const auto t12 = df;
    const auto t13 = htt_r;
    const auto t15 = t12 - epsilon * t13;
    const auto t16 = hrr;
    const auto t20 = static_cast<scalar_type>(1) / (static_cast<scalar_type>(1) + epsilon * t16 * t1);
    const auto t24 = r * r;
    const auto t25 = hee;
    const auto t28 = static_cast<scalar_type>(1) / (t24 + epsilon * t25);
    const auto t30 = htt_theta;
    const auto t35 = t15 * t5 / static_cast<scalar_type>(2);
    const auto t36 = t1 * t20;
    const auto t37 = hrr_t;
    const auto t38 = epsilon * t37;
    const auto t40 = t36 * t38 / static_cast<scalar_type>(2);
    const auto t42 = t6 * t30 / static_cast<scalar_type>(2);
    const auto t43 = t28 * epsilon;
    const auto t44 = hee_t;
    const auto t46 = t43 * t44 / static_cast<scalar_type>(2);
    const auto t47 = std::sin(theta);
    const auto t48 = t47 * t47;
    const auto t50 = hpp;
    const auto t53 = static_cast<scalar_type>(1) / (t24 * t48 + epsilon * t50);
    const auto t55 = hpp_t;
    const auto t57 = t53 * epsilon * t55 / static_cast<scalar_type>(2);
    const auto t59 = t5 / c2;
    const auto t64 = hrr_r;
    const auto t66 = t1 * t1;
    const auto t71 = hrr_theta;
    const auto t76 = t36 * epsilon * t71 / static_cast<scalar_type>(2);
    const auto t78 = hee_r;
    const auto t80 = static_cast<scalar_type>(2) * r + epsilon * t78;
    const auto t82 = t80 * t28 / static_cast<scalar_type>(2);
    const auto t85 = hpp_r;
    const auto t87 = static_cast<scalar_type>(2) * r * t48 + epsilon * t85;
    const auto t89 = t87 * t53 / static_cast<scalar_type>(2);
    const auto t96 = hee_theta;
    const auto t100 = std::cos(theta);
    const auto t103 = hpp_theta;
    const auto t105 = static_cast<scalar_type>(2) * t24 * t47 * t100 + epsilon * t103;
    const auto t107 = t105 * t53 / static_cast<scalar_type>(2);

    christoffel_symbols_type symbols;
    symbols.setZero();
    
    symbols(0, 0, 0) = -t6 * t7 / static_cast<scalar_type>(2);
    symbols(0, 0, 1) = t1 * c2 * t15 * t20 / static_cast<scalar_type>(2);
    symbols(0, 0, 2) = -t28 * c2 * epsilon * t30 / static_cast<scalar_type>(2);
    symbols(0, 1, 0) = t35;
    symbols(0, 1, 1) = t40;
    symbols(0, 2, 0) = -t42;
    symbols(0, 2, 2) = t46;
    symbols(0, 3, 3) = t57;
    
    symbols(1, 0, 0) = t35;
    symbols(1, 0, 1) = t40;
    symbols(1, 1, 0) = t59 * t38 / static_cast<scalar_type>(2);
    symbols(1, 1, 1) = static_cast<scalar_type>(1) / t1 * t20 * (-t12 + epsilon * t64 * t66) / static_cast<scalar_type>(2);
    symbols(1, 1, 2) = -t43 * t71 / static_cast<scalar_type>(2);
    symbols(1, 2, 1) = t76;
    symbols(1, 2, 2) = t82;
    symbols(1, 3, 3) = t89;
    
    symbols(2, 0, 0) = -t42;
    symbols(2, 0, 2) = t46;
    symbols(2, 1, 1) = t76;
    symbols(2, 1, 2) = t82;
    symbols(2, 2, 0) = t59 * epsilon * t44 / static_cast<scalar_type>(2);
    symbols(2, 2, 1) = -t1 * t80 * t20 / static_cast<scalar_type>(2);
    symbols(2, 2, 2) = t43 * t96 / static_cast<scalar_type>(2);
    symbols(2, 3, 3) = t107;
    
    symbols(3, 0, 3) = t57;
    symbols(3, 1, 3) = t89;
    symbols(3, 2, 3) = t107;
    symbols(3, 3, 0) = t59 * epsilon * t55 / static_cast<scalar_type>(2);
    symbols(3, 3, 1) = -t1 * t87 * t20 / static_cast<scalar_type>(2);
    symbols(3, 3, 2) = -t105 * t28 / static_cast<scalar_type>(2);
    
    return symbols;
  }
  
  scalar_type mass = static_cast<scalar_type>(1);
  scalar_type epsilon = static_cast<scalar_type>(0.1);  // Perturbation amplitude
  scalar_type sigma = static_cast<scalar_type>(1);      // Wave frequency
  int l = 0;  // Legendre polynomial order (0 or 1 only)
};
}
