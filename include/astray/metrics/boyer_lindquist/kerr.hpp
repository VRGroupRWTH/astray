#pragma once

#include <cmath>

#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>
#include <astray/math/ipow.hpp>

namespace ast::metrics
{
template <
  typename scalar_type              , 
  typename vector_type              = vector4  <scalar_type>, 
  typename christoffel_symbols_type = tensor444<scalar_type>>
class kerr : public metric<coordinate_system_type::boyer_lindquist, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  using consts = constants<scalar_type>;

  __device__ scalar_type              coordinate_system_parameter() const override
  {
    return angular_momentum / mass;
  }
  __device__ termination_reason       check_termination          (const vector_type& position, const vector_type& direction) const override
  {
    const auto mass_sq = ipow<2>(mass);
    const auto ang_mom_sq = ipow<2>(angular_momentum);
    const auto event_horizon = mass + std::sqrt(mass_sq - ang_mom_sq);
    if (position[1] < static_cast<scalar_type>(0) || position[1] <= (static_cast<scalar_type>(1) + consts::epsilon) * event_horizon)
      return termination_reason::spacetime_breakdown;
    return termination_reason::none;
  }

  __device__ christoffel_symbols_type christoffel_symbols        (const vector_type& position) const override
  {
    const auto r_sq = ipow<2>(position[1]);
    const auto t2   = mass * position[1];
    const auto a_sq = ipow<2>(angular_momentum);
    const auto t5   = r_sq - static_cast<scalar_type>(2) * t2 + a_sq;
    const auto t6   = std::cos(position[2]);
    const auto t6_sq = ipow<2>(t6);
    const auto t8   = a_sq * t6_sq;
    const auto t9   = r_sq + t8;
    const auto t9_sq = ipow<2>(t9);
    const auto t12  = static_cast<scalar_type>(1) / t9_sq / t9;
    const auto t14  = -r_sq + t8;
    const auto t20  = std::sin(position[2]);
    const auto t21  = a_sq * t6 * t20;
    const auto t24  = a_sq + r_sq;
    const auto t26  = static_cast<scalar_type>(1) / t9;
    const auto a_sq_sq = ipow<2>(a_sq);
    const auto t29  = a_sq_sq * t6_sq;
    const auto t30  = r_sq * a_sq;
    const auto t31  = t30 * t6_sq;
    const auto t32  = a_sq * mass;
    const auto t35  = static_cast<scalar_type>(2) * t32 * position[1] * t6_sq;
    const auto r_sq_sq = ipow<2>(r_sq);
    const auto t37  = r_sq * position[1];
    const auto t38  = mass * t37;
    const auto t41  = static_cast<scalar_type>(1) / (t29 + t31 + t30 - t35 + r_sq_sq - static_cast<scalar_type>(2) * t38);
    const auto t42  = t14 * t26 * t41;
    const auto t43  = t24 * mass * t42;
    const auto t44  = mass * angular_momentum;
    const auto t45  = t44 * t42;
    const auto t46  = static_cast<scalar_type>(1) / t9_sq;
    const auto t49  = static_cast<scalar_type>(2) * t2 * t46 * t21;
    const auto t50  = t44 * position[1];
    const auto t51  = static_cast<scalar_type>(1) / t20;
    const auto t55  = static_cast<scalar_type>(2) * t50 * t6 * t51 * t46;
    const auto t58  = -static_cast<scalar_type>(1) + t6_sq;
    const auto t61  = t5 * mass * angular_momentum * t58 * t14 * t12;
    const auto t62  = t20 * t6;
    const auto t66  = static_cast<scalar_type>(2) * t50 * t62 * t24 * t12;
    const auto t68  = static_cast<scalar_type>(1) / t5 * t26;
    const auto t69  = mass * r_sq;
    const auto t77  = t26 * a_sq * t62;
    const auto t78  = t26 * position[1];
    const auto t85  = (t29 - t31 - t30 - static_cast<scalar_type>(3) * r_sq_sq) * mass * angular_momentum * t58 * t26 * t41;
    const auto t6_sq_sq = ipow<2>(t6_sq);
    const auto t88  = position[1] * a_sq_sq * t6_sq_sq;
    const auto t89  = mass * a_sq_sq;
    const auto t90  = t89 * t6_sq;
    const auto t91  = t89 * t6_sq_sq;
    const auto t92  = t69 * t8;
    const auto t95  = static_cast<scalar_type>(2) * t37 * a_sq * t6_sq;
    const auto t96  = t32 * r_sq;
    const auto t97  = r_sq_sq * position[1];
    const auto t102 = (t88 + t90 - t91 - t92 + t95 - t96 + t97 - static_cast<scalar_type>(2) * mass * r_sq_sq) * t26 * t41;
    const auto t112 = static_cast<scalar_type>(2) * t20 * t58 * mass * a_sq * angular_momentum * position[1] * t6 * t46;
    const auto t113 = t6_sq_sq * a_sq_sq;
    const auto t120 = t6 * (t113 + static_cast<scalar_type>(2) * t32 * position[1] - t35 + static_cast<scalar_type>(2) * t31 + r_sq_sq) * t51 * t46;
    const auto t126 = r_sq_sq * a_sq;
    const auto t129 = r_sq * a_sq_sq;

    christoffel_symbols_type symbols;
    symbols.setZero();
    symbols(0, 0, 1) = -t5  * t12 * mass * t14;
    symbols(0, 0, 2) = -static_cast<scalar_type>(2) * t12 * mass * position[1] * t21;
    symbols(0, 1, 0) = -t43;
    symbols(0, 1, 3) = -t45;
    symbols(0, 2, 0) = -t49;
    symbols(0, 2, 3) = -t55;
    symbols(0, 3, 1) = -t61;
    symbols(0, 3, 2) =  t66;
    symbols(1, 0, 0) = -t43;
    symbols(1, 0, 3) = -t45;
    symbols(1, 1, 1) = -t68 * (t69 - position[1] * t4 + t8 * position[1] - t8 * mass);
    symbols(1, 1, 2) =  t68 * t21;
    symbols(1, 2, 1) = -t77;
    symbols(1, 2, 2) =  t78;
    symbols(1, 3, 0) = -t85;
    symbols(1, 3, 3) = t102;
    symbols(2, 0, 0) = -t49;
    symbols(2, 0, 3) = -t55;
    symbols(2, 1, 1) = -t77;
    symbols(2, 1, 2) =  t78;
    symbols(2, 2, 1) = -t5 * t26 * position[1];
    symbols(2, 2, 2) = -t77;
    symbols(2, 3, 0) = -t112;
    symbols(2, 3, 3) =  t120;
    symbols(3, 0, 1) = -t61;
    symbols(3, 0, 2) =  t66;
    symbols(3, 1, 0) = -t85;
    symbols(3, 1, 3) =  t102;
    symbols(3, 2, 0) = -t112;
    symbols(3, 2, 3) =  t120;
    symbols(3, 3, 1) =  t5 * t58 * (t97 + t95 + t88 - t96 + t92 + t90 - t91) * t12;
    symbols(3, 3, 2) = -t62 * (t36 * t1 + static_cast<scalar_type>(2) * t126 * t7 + t129 * t87 + t126
      + static_cast<scalar_type>(2) * t129 * t7   + t28 * t4 * t87 
      + static_cast<scalar_type>(4) * t38  * t4   - static_cast<scalar_type>(4) * t38 * t8 
      - static_cast<scalar_type>(2) * t2   * t113 + static_cast<scalar_type>(2) * t89 * position[1]) * t12;
    return symbols;
  }
  
  scalar_type mass             = static_cast<scalar_type>(1);
  scalar_type angular_momentum = static_cast<scalar_type>(1);
};
}