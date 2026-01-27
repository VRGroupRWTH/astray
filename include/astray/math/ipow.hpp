#pragma once

namespace ast
{
// Compile-time integer power function for improved code clarity and intent.
// Computes base^N where N is a compile-time constant integer.
// More efficient than std::pow for small integer powers while maintaining intent.
template <int N, typename T>
__device__ __host__ constexpr T ipow(const T& base)
{
  static_assert(N >= 0, "ipow only supports non-negative integer powers");
  
  if constexpr (N == 0)
    return T(1);
  else if constexpr (N == 1)
    return base;
  else if constexpr (N == 2)
    return base * base;
  else if constexpr (N == 3)
    return base * base * base;
  else if constexpr (N == 4)
  {
    const T base_sq = base * base;
    return base_sq * base_sq;
  }
  else if constexpr (N % 2 == 0)
  {
    const T half_power = ipow<N / 2>(base);
    return half_power * half_power;
  }
  else
  {
    return base * ipow<N - 1>(base);
  }
}
}
