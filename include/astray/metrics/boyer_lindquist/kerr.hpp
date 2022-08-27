#pragma once

#include <cmath>

#include <astray/core/metric.hpp>
#include <astray/math/constants.hpp>

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
    const auto event_horizon = mass + std::sqrt(static_cast<scalar_type>(std::pow(mass, 2)) - static_cast<scalar_type>(std::pow(angular_momentum, 2)));
    if (position[1] < static_cast<scalar_type>(0) || position[1] <= (static_cast<scalar_type>(1) + consts::epsilon) * event_horizon)
      return termination_reason::spacetime_breakdown;
    return termination_reason::none;
  }

  __device__ christoffel_symbols_type christoffel_symbols        (const vector_type& position) const override
  {
    const auto t1   = static_cast<scalar_type>(std::pow(position[1], 2));
    const auto t2   = mass * position[1];
    const auto t4   = static_cast<scalar_type>(std::pow(angular_momentum, 2));
    const auto t5   = t1 - static_cast<scalar_type>(2) * t2 + t4;
    const auto t6   = std::cos(position[2]);
    const auto t7   = static_cast<scalar_type>(std::pow(t6, 2));
    const auto t8   = t4 * t7;
    const auto t9   = t1 + t8;
    const auto t10  = static_cast<scalar_type>(std::pow(t9, 2));
    const auto t12  = static_cast<scalar_type>(1) / t10 / t9;
    const auto t14  = -t1 + t8;
    const auto t20  = std::sin(position[2]);
    const auto t21  = t4 * t6 * t20;
    const auto t24  = t4 + t1;
    const auto t26  = static_cast<scalar_type>(1) / t9;
    const auto t28  = static_cast<scalar_type>(std::pow(t4, 2));
    const auto t29  = t28 * t7;
    const auto t30  = t1 * t4;
    const auto t31  = t30 * t7;
    const auto t32  = t4 * mass;
    const auto t35  = static_cast<scalar_type>(2) * t32 * position[1] * t7;
    const auto t36  = static_cast<scalar_type>(std::pow(t1, 2));
    const auto t37  = t1 * position[1];
    const auto t38  = mass * t37;
    const auto t41  = static_cast<scalar_type>(1) / (t29 + t31 + t30 - t35 + t36 - static_cast<scalar_type>(2) * t38);
    const auto t42  = t14 * t26 * t41;
    const auto t43  = t24 * mass * t42;
    const auto t44  = mass * angular_momentum;
    const auto t45  = t44 * t42;
    const auto t46  = static_cast<scalar_type>(1) / t10;
    const auto t49  = static_cast<scalar_type>(2) * t2 * t46 * t21;
    const auto t50  = t44 * position[1];
    const auto t51  = static_cast<scalar_type>(1) / t20;
    const auto t55  = static_cast<scalar_type>(2) * t50 * t6 * t51 * t46;
    const auto t58  = -static_cast<scalar_type>(1) + t7;
    const auto t61  = t5 * mass * angular_momentum * t58 * t14 * t12;
    const auto t62  = t20 * t6;
    const auto t66  = static_cast<scalar_type>(2) * t50 * t62 * t24 * t12;
    const auto t68  = static_cast<scalar_type>(1) / t5 * t26;
    const auto t69  = mass * t1;
    const auto t77  = t26 * t4 * t62;
    const auto t78  = t26 * position[1];
    const auto t85  = (t29 - t31 - t30 - static_cast<scalar_type>(3) * t36) * mass * angular_momentum * t58 * t26 * t41;
    const auto t87  = static_cast<scalar_type>(std::pow(t7, 2));
    const auto t88  = position[1] * t28 * t87;
    const auto t89  = mass * t28;
    const auto t90  = t89 * t7;
    const auto t91  = t89 * t87;
    const auto t92  = t69 * t8;
    const auto t95  = static_cast<scalar_type>(2) * t37 * t4 * t7;
    const auto t96  = t32 * t1;
    const auto t97  = t36 * position[1];
    const auto t102 = (t88 + t90 - t91 - t92 + t95 - t96 + t97 - static_cast<scalar_type>(2) * mass * t36) * t26 * t41;
    const auto t112 = static_cast<scalar_type>(2) * t20 * t58 * mass * t4 * angular_momentum * position[1] * t6 * t46;
    const auto t113 = t87 * t28;
    const auto t120 = t6 * (t113 + static_cast<scalar_type>(2) * t32 * position[1] - t35 + static_cast<scalar_type>(2) * t31 + t36) * t51 * t46;
    const auto t126 = t36 * t4;
    const auto t129 = t1 * t28;

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