#pragma once

#include <cstddef>
#include <cstdint>
#include <iterator>

#include <astray/math/ode/error/controller/proportional_integral_controller.hpp>
#include <astray/math/ode/tableau/tableau_traits.hpp>

namespace ast
{
template <
  typename method_type_          , 
  typename problem_type_         , 
  typename error_evaluator_type_ = proportional_integral_controller<std::remove_reference_t<typename problem_type_::time_type>, typename method_type_::tableau_type>>
class adaptive_step_iterator
{
public:
  using method_type          = method_type_;
  using problem_type         = problem_type_;
  using error_evaluator_type = error_evaluator_type_;
  using time_type            = typename problem_type::time_type;

  using iterator_category    = std::input_iterator_tag; // Single pass read forward.
  using difference_type      = std::ptrdiff_t;
  using value_type           = problem_type_;
  using pointer              = value_type const*;
  using reference            = value_type const&;

  __device__ constexpr reference               operator* () const noexcept
  {
    return  problem;
  }
  __device__ constexpr pointer                 operator->() const noexcept
  {
    return &problem;
  }

  __device__ constexpr adaptive_step_iterator& operator++()
  {
    if constexpr (is_extended_butcher_tableau_v<typename method_type::tableau_type>)
    {
      const auto result     = method_type::apply      (problem, step_size);
      const auto evaluation = error_evaluator.evaluate(problem, step_size, result);

      if (evaluation.accept)
      {
        problem.value = result.value;
        problem.time += step_size;
        step_size     = evaluation.next_step_size;
      }
      else
      {
        step_size     = evaluation.next_step_size;
        return operator++(); // Retry with the adapted step size.
      }
    }
    else
    {
      // Non-extended Butcher tableau may be used with adaptive step iterators, but their step size will not be adapted.
      problem.value = method_type::apply(problem, step_size);
      problem.time += step_size;
    }
    return *this;
  }
  __device__ constexpr adaptive_step_iterator  operator++(std::int32_t)
  {
    adaptive_step_iterator temp = *this;
    ++(*this);
    return temp;
  }

  __device__ friend constexpr bool             operator==(const adaptive_step_iterator& lhs, const adaptive_step_iterator& rhs) noexcept
  {
    return lhs.problem == rhs.problem && lhs.step_size == rhs.step_size && lhs.error_evaluator == rhs.error_evaluator;
  }
  
  problem_type         problem        ;
  time_type            step_size      ;
  error_evaluator_type error_evaluator;
};
}