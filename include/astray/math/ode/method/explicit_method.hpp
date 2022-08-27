#pragma once

#include <array>
#include <utility>

#include <astray/math/ode/error/extended_result.hpp>
#include <astray/math/ode/tableau/tableau_traits.hpp>
#include <astray/math/ode/utility/constexpr_for.hpp>
#include <astray/math/ode/utility/triangular_number.hpp>

namespace ast
{
template <typename tableau_type_>
class explicit_method
{
public:
  using tableau_type = tableau_type_;

  template <typename problem_type>
  __device__ static constexpr auto apply(const problem_type& problem, const typename problem_type::time_type step_size)
  {
    using value_type = typename problem_type::function_type::result_type;

    std::array<value_type, stages_v<tableau_type>> stages;
    constexpr_for<0, stages_v<tableau_type>, 1>([&problem, &step_size, &stages] (auto i)
    {
      value_type sum;
      constexpr_for<0, i.value, 1>([&stages, &sum, &i] (auto j)
      {
        sum += std::get<j>(stages) * std::get<triangular_number<i.value - 1> + j.value>(tableau_a<tableau_type>);
      });
      std::get<i>(stages) = problem.function(problem.time + std::get<i>(tableau_c<tableau_type>) * step_size, problem.value + sum * step_size);
    });

    if constexpr (is_extended_butcher_tableau_v<tableau_type>)
    {
      value_type higher, lower;
      constexpr_for<0, stages_v<tableau_type>, 1>([&stages, &higher, &lower] (auto i)
      {
        higher += std::get<i>(stages) * std::get<i>(tableau_b <tableau_type>);
        lower  += std::get<i>(stages) * std::get<i>(tableau_bs<tableau_type>);
      });
      return extended_result<value_type> {problem.value + higher * step_size, (higher - lower) * step_size};
    }
    else
    {
      value_type sum;
      constexpr_for<0, stages_v<tableau_type>, 1>([&stages, &sum] (auto i)
      {
        sum += std::get<i>(stages) * std::get<i>(tableau_b<tableau_type>);
      });
      return value_type(problem.value + sum * step_size);
    }
  }
};
}