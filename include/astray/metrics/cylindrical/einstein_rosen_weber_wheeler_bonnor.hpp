#pragma once

#include <cmath>

#include <astray/core/metric.hpp>

namespace ast::metrics
{
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class einstein_rosen_weber_wheeler_bonnor : public metric<coordinate_system_type::cylindrical, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    scalar_type psi  , psi_dt  , psi_dr  ;
    scalar_type gamma, gamma_dt, gamma_dr;

    christoffel_symbols_type symbols;
    symbols.setZero();

    {
      const auto t1  = std::sqrt(static_cast<scalar_type>(2));
      const auto t2  = t1 * constant;
      const auto t3  = static_cast<scalar_type>(std::pow(pulse_width, 2));
      const auto t4  = static_cast<scalar_type>(std::pow(t3, 2));
      const auto t5  = static_cast<scalar_type>(std::pow(position[1], 2));
      const auto t6  = t3 * t5;
      const auto t8  = static_cast<scalar_type>(std::pow(position[0], 2));
      const auto t9  = t3 * t8;
      const auto t11 = static_cast<scalar_type>(std::pow(t5, 2));
      const auto t14 = static_cast<scalar_type>(std::pow(t8, 2));
      const auto t15 = t4 + static_cast<scalar_type>(2) * t6 + static_cast<scalar_type>(2) * t9 + t11 - static_cast<scalar_type>(2) * t5 * t8 + t14;
      const auto t16 = std::sqrt(t15);
      const auto t17 = t16 + t3 + t5 - t8;
      const auto t18 = t3 + t5 - t8;
      const auto t19 = static_cast<scalar_type>(std::pow(t18, 2));
      const auto t20 = static_cast<scalar_type>(4) * t9;
      const auto t21 = t19 + t20;
      const auto t22 = static_cast<scalar_type>(1) / t21;
      const auto t24 = std::sqrt(t17 * t22);
      const auto t25 = static_cast<scalar_type>(1) / t24;
      const auto t26 = static_cast<scalar_type>(1) / t16;
      const auto t27 = t3 * position[0];
      const auto t30 = t27 - t5 * position[0] + t8 * position[0];
      const auto t36 = static_cast<scalar_type>(std::pow(t21, 2));
      const auto t37 = static_cast<scalar_type>(1) / t36;
      const auto t38 = t17 * t37;
      const auto t40 = static_cast<scalar_type>(4) * t18 * position[0];
      const auto t41 = static_cast<scalar_type>(8) * t27;
      const auto t42 = -t40 + t41;
      const auto t48 = t3 * position[1];
      const auto t49 = t5 * position[1];
      const auto t51 = t48 + t49 - position[1] * t8;
      const auto t64 = static_cast<scalar_type>(std::pow(constant, 2));
      const auto t66 = t64 / t3;
      const auto t71 = t19 - t20;
      const auto t74 = t71 / t36 / t21;
      const auto t83 = (t5 - t3 - t8) / t16 / t15;
      const auto t92 = t3 * t49;
      
      psi      = t1 * constant * t24;
      gamma    = t64 / t3 * (static_cast<scalar_type>(1) - static_cast<scalar_type>(2) * t6 * (t19 - t20) / t36 + (t5 - t3 - t8) / t16) / static_cast<scalar_type>(2);

      psi_dt   = t2 * t25 * ((static_cast<scalar_type>(2) * t26 * t30 - static_cast<scalar_type>(2) * position[0]) * 
        t22 - t38 * t42) / static_cast<scalar_type>(2);
      psi_dr   = t2 * t25 * ((static_cast<scalar_type>(2) * t26 * t51 + static_cast<scalar_type>(2) * position[1]) * 
        t22 - static_cast<scalar_type>(4) * t38 * t18 * position[1]) / static_cast<scalar_type>(2);

      gamma_dt = t66 * (-static_cast<scalar_type>(2) * t6 * (-t40 - t41) * t37 + static_cast<scalar_type>(4) * t6 * t74 * 
        t42 - static_cast<scalar_type>(2) * position[0] * t26 - static_cast<scalar_type>(2) * t83 * t30) / static_cast<scalar_type>(2);
      gamma_dr = t66 * (-static_cast<scalar_type>(4) * t48 * t71 * t37 - static_cast<scalar_type>(8) * t92 * t18 * t37 + static_cast<scalar_type>(16) * 
        t92 * t74 * t18 + static_cast<scalar_type>(2) * position[1] * t26 - static_cast<scalar_type>(2) * t83 * t51) / static_cast<scalar_type>(2);
    }

    {
      const auto t3  = gamma_dt - psi_dt;
      const auto t6  = gamma_dr - psi_dr;
      const auto t9  = static_cast<scalar_type>(-1) + position[1] * psi_dr;
      const auto t10 = t9 / position[1];
      const auto t12 = std::exp(gamma);
      const auto t13 = static_cast<scalar_type>(std::pow(t12, 2));
      const auto t14 = static_cast<scalar_type>(1) / t13;
      const auto t15 = static_cast<scalar_type>(std::pow(position[1], 2));
      const auto t20 = psi;
      const auto t21 = std::exp(t20);
      const auto t22 = static_cast<scalar_type>(std::pow(t21, 2));
      const auto t23 = static_cast<scalar_type>(std::pow(t22, 2));
      const auto t24 = t14 * t23;

      symbols(0, 0, 0) =  t3;
      symbols(0, 0, 1) =  t6;
      symbols(0, 1, 0) =  t6;
      symbols(0, 1, 1) =  t3;
      symbols(0, 2, 2) = -psi_dt;
      symbols(0, 3, 3) =  psi_dt;
      symbols(1, 0, 0) =  t6;
      symbols(1, 0, 1) =  t3;
      symbols(1, 1, 0) =  t3;
      symbols(1, 1, 1) =  t6;
      symbols(1, 2, 2) = -t10;
      symbols(1, 3, 3) =  psi_dr;
      symbols(2, 0, 2) = -psi_dt;
      symbols(2, 1, 2) = -t10;
      symbols(2, 2, 0) = -t14 * t15 * psi_dt;
      symbols(2, 2, 1) =  t14 * position[1] * t9;
      symbols(3, 0, 3) =  psi_dt;
      symbols(3, 1, 3) =  psi_dr;
      symbols(3, 3, 0) =  t24 * psi_dt;
      symbols(3, 3, 1) = -t24 * psi_dr;
    }

    return symbols;
  }

  scalar_type pulse_width = static_cast<scalar_type>(1);
  scalar_type constant    = static_cast<scalar_type>(1);
};
}