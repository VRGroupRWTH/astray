#pragma once

#include <cmath>

#include <astray/core/metric.hpp>
#include <astray/parallel/thrust.hpp>

#if (THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA)
#include <cuda_runtime_api.h>
#endif

namespace ast::metrics
{
template <
  typename scalar_type              , 
  typename vector_type              = ast::vector4  <scalar_type>, 
  typename christoffel_symbols_type = ast::tensor444<scalar_type>>
class bessel : public metric<coordinate_system_type::cartesian, scalar_type, vector_type, christoffel_symbols_type>
{
public:
  __device__ christoffel_symbols_type christoffel_symbols(const vector_type& position) const override
  {
    const auto t1   = static_cast<scalar_type>(std::pow(position[1], 2));
    const auto t2   = static_cast<scalar_type>(std::pow(position[2], 2));
    const auto t3   = t1 + t2;
    const auto t4   = std::sqrt(t3);
    const auto t5   = bessel::cyl_bessel_j0(t4);
    const auto t6   = C * t5;
    const auto t7   = std::sin(position[0]);
    const auto t8   = C * t4;
    const auto t9   = bessel::cyl_bessel_j1(t4);
    const auto t10  = std::cos(position[0]);
    const auto t11  = t9 * t10;
    const auto t18  = static_cast<scalar_type>(std::pow(t5 , 2));
    const auto t20  = static_cast<scalar_type>(std::pow(t9 , 2));
    const auto t21  = static_cast<scalar_type>(std::pow(t10, 2));
    const auto t27  = static_cast<scalar_type>(1) / t4;
    const auto t28  = (t8 * t18 + t8 * t20 * t21 - t8 * t18 * t21 + t11) * t27;
    const auto t29  = C * position[1] * t28;
    const auto t30  = C * position[2];
    const auto t31  = t28 * t30;
    const auto t32  = C * t1;
    const auto t37  = t6 * t27;
    const auto t38  = (t4 + static_cast<scalar_type>(2) * t32 * t11) * t7 * t37;
    const auto t39  = t7 * t5;
    const auto t41  = static_cast<scalar_type>(std::pow(C, 2));
    const auto t43  = position[1] * t27;
    const auto t46  = static_cast<scalar_type>(2) * t11 * t39 * t41 * position[2] * t43;
    const auto t47  = C * t9;
    const auto t53  = (static_cast<scalar_type>(2) * t47 * t10 * t2 + t4) * t7 * t37;
    const auto t54  = t6 * t7;
    const auto t64  = t41 * (t18 * t1 + t18 * t2 + t20 * t1 + t20 * t2 - static_cast<scalar_type>(2) * t4 * t5 * t9 * t21);
    const auto t65  = std::exp(-t64);
    const auto t67  = C * t65 * t5;
    const auto t68  = std::exp(t64);
    const auto t69  = t1 * t1;
    const auto t70  = t68 * t69;
    const auto t71  = t47 * t10;
    const auto t72  = t70 * t71;
    const auto t73  = static_cast<scalar_type>(2) * t72;
    const auto t74  = t68 * t1;
    const auto t76  = t11 * t2;
    const auto t77  = C * t74 * t76;
    const auto t78  = static_cast<scalar_type>(2) * t77;
    const auto t79  = t74 * t4;
    const auto t80  = t2 * t4;
    const auto t84  = static_cast<scalar_type>(1) / t4 / t3;
    const auto t87  = t2 * t2;
    const auto t88  = t87 * t68;
    const auto t89  = t88 * t71;
    const auto t90  = static_cast<scalar_type>(2) * t89;
    const auto t92  = t47 * t10 * t87;
    const auto t93  = t68 * t4;
    const auto t94  = t93 * t2;
    const auto t95  = t74 * t41;
    const auto t96  = t4 * t20;
    const auto t97  = t21 * t2;
    const auto t99  = t95 * t96 * t97;
    const auto t100 = t32 * t76;
    const auto t101 = t4 * t18;
    const auto t103 = t95 * t101 * t2;
    const auto t105 = t95 * t101 * t97;
    const auto t106 = static_cast<scalar_type>(3) * t77;
    const auto t108 = t41 * t4 * t18;
    const auto t109 = t70 * t108;
    const auto t110 = t70 * t41;
    const auto t111 = t96 * t21;
    const auto t112 = t110 * t111;
    const auto t113 = t101 * t21;
    const auto t114 = t110 * t113;
    const auto t115 = t90 - t92 + t94 + t99 - t80 - t100 + t103 - t105 + t106 + t72 + t109 + t112 - t114;
    const auto t117 = t3 * t3;
    const auto t119 = static_cast<scalar_type>(1) / t4 / t117;
    const auto t120 = position[1] * t119;
    const auto t122 = t92 - t99 - t94 + t100 + t105 - t103 + t80 + t77 + t72 - t109 + t114 - t112;
    const auto t124 = position[2] * t119;
    const auto t136 = t65 * position[1] * t30 * t39 * (static_cast<scalar_type>(2) * t74 * t71 + static_cast<scalar_type>(2) * t68 * C * t76 + t93 - t4) * t84;
    const auto t138 = t47 * t10 * t69;
    const auto t139 = t4 * t1;
    const auto t140 = t89 + t99 + t100 - t105 + t77 + t103 + t138 - t114 + t139 + t112 + t109 - t79;
    const auto t142 = t140 * t65 * t124;
    const auto t143 = t88 * t108;
    const auto t144 = t88 * t41;
    const auto t145 = t144 * t113;
    const auto t146 = t144 * t111;
    const auto t147 = t92 + t143 - t145 + t146 - t105 + t99 + t80 + t100 + t103 + t77 - t94 + t72;
    const auto t149 = t147 * t65 * t120;
    const auto t150 = t43 * t10;
    const auto t151 = t47 * t150;
    const auto t156 = t145 - t143 - t146 + t89 - t103 - t99 + t77 + t105 + t100 + t139 + t138 - t79;
    const auto t159 = -t89 - t143 - t146 + t145 - t99 + t105 + t100 - t103 - t106 + t138 - t73 - t79 + t139;
    const auto t164 = t47 * t27 * position[2] * t10;
    const auto t165 = C * t18;
    const auto t168 = C * t20;
    const auto t179 = std::exp(-C * (t165 * t1 + t165 * t2 + t168 * t1 + t168 * t2 - static_cast<scalar_type>(2) * t8 * t5 * t9 * t21 - static_cast<scalar_type>(4) * t5 * t10));
    const auto t180 = C * t179;

    christoffel_symbols_type symbols;
    symbols.setZero();
    symbols(0, 0, 0) =  t6 * t7 * (static_cast<scalar_type>(2) * t8 * t11 + static_cast<scalar_type>(1));
    symbols(0, 0, 1) =  t29;
    symbols(0, 0, 2) =  t31;
    symbols(0, 1, 0) =  t29;
    symbols(0, 1, 1) =  t38;
    symbols(0, 1, 2) =  t46;
    symbols(0, 2, 0) =  t31;
    symbols(0, 2, 1) =  t46;
    symbols(0, 2, 2) =  t53;
    symbols(0, 3, 3) = -t54;
    symbols(1, 0, 0) =  t29;
    symbols(1, 0, 1) =  t38;
    symbols(1, 0, 2) =  t46;
    symbols(1, 1, 0) =  t67 * t7 * (t73 + t78 + t79 + t80) * t84;
    symbols(1, 1, 1) =  t115 * t65 * t120;
    symbols(1, 1, 2) = -t122 * t65 * t124;
    symbols(1, 2, 0) =  t136;
    symbols(1, 2, 1) =  t142;
    symbols(1, 2, 2) =  t149;
    symbols(1, 3, 3) = -t151;
    symbols(2, 0, 0) =  t31;
    symbols(2, 0, 1) =  t46;
    symbols(2, 0, 2) =  t53;
    symbols(2, 1, 0) =  t136;
    symbols(2, 1, 1) =  t142;
    symbols(2, 1, 2) =  t149;
    symbols(2, 2, 0) =  t67 * t7 * (t78 + t90 + t94 + t139) * t84;
    symbols(2, 2, 1) = -t156 * t65 * t120;
    symbols(2, 2, 2) = -t159 * t65 * t124;
    symbols(2, 3, 3) = -t164;
    symbols(3, 0, 3) = -t54;
    symbols(3, 1, 3) = -t151;
    symbols(3, 2, 3) = -t164;
    symbols(3, 3, 0) = -t180 * t39;
    symbols(3, 3, 1) =  t180 * t9 * t150;
    symbols(3, 3, 2) =  t179 * position[2] * C * t11 * t27;
    return symbols;
  }

  __device__ static scalar_type       cyl_bessel_j0      (const scalar_type x)
  {
#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
      return static_cast<scalar_type>(j0(x));
#else
      return std::cyl_bessel_j(0, x);
#endif
  }
  __device__ static scalar_type       cyl_bessel_j1      (const scalar_type x)
  {
#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
      return static_cast<scalar_type>(j1(x));
#else
      return std::cyl_bessel_j(1, x);
#endif
  }

  scalar_type C = static_cast<scalar_type>(1);
};
}