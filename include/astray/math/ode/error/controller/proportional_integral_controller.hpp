#pragma once

#include <algorithm>
#include <cmath>
#include <complex>

#include <astray/math/ode/algebra/quantity_operations.hpp>
#include <astray/math/ode/error/error_evaluation.hpp>
#include <astray/math/ode/tableau/tableau_traits.hpp>

namespace ast
{
template <typename type, typename tableau_type>
struct proportional_integral_controller
{
  // Reference: https://doi.org/10.1007/978-3-642-05221-7 Chapter: IV.2, Section: A PI Step Size Control, Equation: 2.43c
  template <typename problem_type, typename extended_result_type>
  __device__ constexpr error_evaluation<type> evaluate(const problem_type& problem, const type step_size, const extended_result_type& result)
  {
    using operations = quantity_operations<typename extended_result_type::type>;

    type squared_sum(0);
    operations::for_each([&] (const auto& p, const auto& r, const auto& e)
    {
      squared_sum += static_cast<type>(std::pow(std::abs(e) / (absolute_tolerance + relative_tolerance * std::max(std::abs(p), std::abs(r))), 2));
    }, problem.value, result.value, result.error);

    type error   = std::sqrt(squared_sum / operations::size(problem.value)); // std::real(squared_sum) unavailable in CUDA until C++20 support.
    type optimal = factor * std::pow(type(1) / std::abs(error), alpha) * std::pow(std::abs(previous_error), beta);
    type limited = std::min (factor_maximum, std::max(factor_minimum, optimal));

    const bool accept = error <= type(1);
    if (accept)
      previous_error = error;

    return {accept, step_size * limited};
  }
  
  const type absolute_tolerance = type(1e-6);
  const type relative_tolerance = type(1e-3);
  const type factor             = type(0.8 );
  const type factor_minimum     = type(1e-2);
  const type factor_maximum     = type(1e+2);
  const type alpha              = type(7.0 / (10.0 * order_v<tableau_type>));
  const type beta               = type(4.0 / (10.0 * order_v<tableau_type>));
  type       previous_error     = type(1e-3);
};
}