#pragma once

#include <cmath>

namespace ast::math
{
// CUDA-compatible special function implementations
// These are device-callable approximations suitable for use in metric calculations

// Lambert W function (principal branch W0) using Halley's method
// Converges for x >= -1/e ≈ -0.367879
template <typename scalar_type>
__device__ __host__ inline scalar_type lambert_w0(scalar_type x)
{
  // Handle special cases
  if (x < static_cast<scalar_type>(-0.367879))
    return static_cast<scalar_type>(0); // Outside domain
  if (x == static_cast<scalar_type>(0))
    return static_cast<scalar_type>(0);
  
  // Initial guess using asymptotic approximation
  scalar_type w;
  if (x < static_cast<scalar_type>(1))
  {
    // For small x, use series expansion: W(x) ≈ x - x² + 3x³/2 - 8x⁴/3 + ...
    w = x;
  }
  else
  {
    // For large x, use W(x) ≈ log(x) - log(log(x))
    const auto lx = std::log(x);
    w = lx - std::log(lx);
  }
  
  // Halley's method iteration: w_new = w - (w*e^w - x) / (e^w * (w+1) - (w+2)*(w*e^w - x)/(2*(w+1)))
  // Simplifies to: w_new = w - f/f' where f = w*e^w - x and f' adjusted for Halley
  constexpr int max_iter = 20;
  constexpr scalar_type tol = static_cast<scalar_type>(1e-10);
  
  for (int i = 0; i < max_iter; ++i)
  {
    const auto ew = std::exp(w);
    const auto wew = w * ew;
    const auto f = wew - x;
    
    if (std::abs(f) < tol)
      break;
    
    // Halley's method update
    const auto fp = ew * (w + static_cast<scalar_type>(1));
    const auto fpp = ew * (w + static_cast<scalar_type>(2));
    const auto delta = f / (fp - f * fpp / (static_cast<scalar_type>(2) * fp));
    
    w = w - delta;
  }
  
  return w;
}

// Bessel function J0(x) using polynomial approximation
// Accurate to ~7 decimal places for |x| < 8
template <typename scalar_type>
__device__ __host__ inline scalar_type bessel_j0(scalar_type x)
{
  const auto ax = std::abs(x);
  
  if (ax < static_cast<scalar_type>(8))
  {
    // Polynomial approximation for |x| < 8
    const auto y = x * x;
    const auto ans1 = static_cast<scalar_type>(57568490574.0) +
                      y * (static_cast<scalar_type>(-13362590354.0) +
                      y * (static_cast<scalar_type>(651619640.7) +
                      y * (static_cast<scalar_type>(-11214424.18) +
                      y * (static_cast<scalar_type>(77392.33017) +
                      y * static_cast<scalar_type>(-184.9052456)))));
    const auto ans2 = static_cast<scalar_type>(57568490411.0) +
                      y * (static_cast<scalar_type>(1029532985.0) +
                      y * (static_cast<scalar_type>(9494680.718) +
                      y * (static_cast<scalar_type>(59272.64853) +
                      y * (static_cast<scalar_type>(267.8532712) +
                      y * static_cast<scalar_type>(1.0)))));
    return ans1 / ans2;
  }
  else
  {
    // Asymptotic approximation for |x| >= 8
    const auto z = static_cast<scalar_type>(8) / ax;
    const auto y = z * z;
    const auto xx = ax - static_cast<scalar_type>(0.785398164);
    
    const auto ans1 = static_cast<scalar_type>(1.0) +
                      y * (static_cast<scalar_type>(-0.1098628627e-2) +
                      y * (static_cast<scalar_type>(0.2734510407e-4) +
                      y * (static_cast<scalar_type>(-0.2073370639e-5) +
                      y * static_cast<scalar_type>(0.2093887211e-6))));
    const auto ans2 = static_cast<scalar_type>(-0.1562499995e-1) +
                      y * (static_cast<scalar_type>(0.1430488765e-3) +
                      y * (static_cast<scalar_type>(-0.6911147651e-5) +
                      y * (static_cast<scalar_type>(0.7621095161e-6) +
                      y * static_cast<scalar_type>(-0.934945152e-7))));
    
    return std::sqrt(static_cast<scalar_type>(0.636619772) / ax) *
           (std::cos(xx) * ans1 - z * std::sin(xx) * ans2);
  }
}

// Bessel function J1(x) using polynomial approximation
// Accurate to ~7 decimal places for |x| < 8
template <typename scalar_type>
__device__ __host__ inline scalar_type bessel_j1(scalar_type x)
{
  const auto ax = std::abs(x);
  
  if (ax < static_cast<scalar_type>(8))
  {
    // Polynomial approximation for |x| < 8
    const auto y = x * x;
    const auto ans1 = x * (static_cast<scalar_type>(72362614232.0) +
                      y * (static_cast<scalar_type>(-7895059235.0) +
                      y * (static_cast<scalar_type>(242396853.1) +
                      y * (static_cast<scalar_type>(-2972611.439) +
                      y * (static_cast<scalar_type>(15704.48260) +
                      y * static_cast<scalar_type>(-30.16036606))))));
    const auto ans2 = static_cast<scalar_type>(144725228442.0) +
                      y * (static_cast<scalar_type>(2300535178.0) +
                      y * (static_cast<scalar_type>(18583304.74) +
                      y * (static_cast<scalar_type>(99447.43394) +
                      y * (static_cast<scalar_type>(376.9991397) +
                      y * static_cast<scalar_type>(1.0)))));
    return ans1 / ans2;
  }
  else
  {
    // Asymptotic approximation for |x| >= 8
    const auto z = static_cast<scalar_type>(8) / ax;
    const auto y = z * z;
    const auto xx = ax - static_cast<scalar_type>(2.356194491);
    
    const auto ans1 = static_cast<scalar_type>(1.0) +
                      y * (static_cast<scalar_type>(0.183105e-2) +
                      y * (static_cast<scalar_type>(-0.3516396496e-4) +
                      y * (static_cast<scalar_type>(0.2457520174e-5) +
                      y * static_cast<scalar_type>(-0.240337019e-6))));
    const auto ans2 = static_cast<scalar_type>(0.04687499995) +
                      y * (static_cast<scalar_type>(-0.2002690873e-3) +
                      y * (static_cast<scalar_type>(0.8449199096e-5) +
                      y * (static_cast<scalar_type>(-0.88228987e-6) +
                      y * static_cast<scalar_type>(0.105787412e-6))));
    
    const auto ans = std::sqrt(static_cast<scalar_type>(0.636619772) / ax) *
                     (std::cos(xx) * ans1 - z * std::sin(xx) * ans2);
    return x < static_cast<scalar_type>(0) ? -ans : ans;
  }
}

} // namespace ast::math
