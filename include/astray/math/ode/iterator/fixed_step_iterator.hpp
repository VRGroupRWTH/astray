#pragma once

#include <cstddef>
#include <cstdint>
#include <iterator>

#include <astray/math/ode/tableau/tableau_traits.hpp>

namespace ast
{
template <typename method_type_, typename problem_type_>
class fixed_step_iterator
{
public:
  using method_type       = method_type_;
  using problem_type      = problem_type_;
  using time_type         = typename problem_type::time_type;

  using iterator_category = std::input_iterator_tag; // Single pass read forward.
  using difference_type   = std::ptrdiff_t;
  using value_type        = problem_type_;
  using pointer           = value_type const*;
  using reference         = value_type const&;

  __device__ constexpr reference            operator* () const noexcept
  {
    return problem;
  }
  __device__ constexpr pointer              operator->() const noexcept
  {
    return &problem;
  }

  __device__ constexpr fixed_step_iterator& operator++()
  {
    if constexpr (is_extended_butcher_tableau_v<typename method_type::tableau_type>)
      problem.value = method_type::apply(problem, step_size).value;
    else
      problem.value = method_type::apply(problem, step_size);

    problem.time += step_size;

    return *this;
  }
  __device__ constexpr fixed_step_iterator  operator++(std::int32_t)
  {
    fixed_step_iterator temp = *this;
    ++(*this);
    return temp;
  }

  __device__ friend constexpr bool          operator==(const fixed_step_iterator& lhs, const fixed_step_iterator& rhs) noexcept
  {
    return lhs.problem == rhs.problem && lhs.step_size == rhs.step_size;
  }
  
  problem_type problem  ;
  time_type    step_size;
};
}