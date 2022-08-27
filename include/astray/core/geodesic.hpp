#pragma once

#include <astray/core/termination_reason.hpp>
#include <astray/math/ode/ode.hpp>
#include <astray/math/linear_algebra.hpp>
#include <astray/parallel/thrust.hpp>

namespace ast
{
template <typename scalar_type_, typename tableau_type_, typename error_evaluator_type_ = proportional_integral_controller<scalar_type_, tableau_type_>>
class geodesic
{
public:
  using scalar_type          = scalar_type_;
  using tableau_type         = tableau_type_;
  using bounds_type          = aabb4<scalar_type>;
  using error_evaluator_type = error_evaluator_type_;

  template <typename ray_type, typename metric_type>
  __device__ static constexpr termination_reason integrate(
    ray_type&                   ray             ,
    const metric_type&          metric          ,
    const std::size_t           iterations      , 
    const scalar_type           lambda_step_size, 
    const scalar_type           lambda          = static_cast<scalar_type>(0),
    const bounds_type&          bounds          = bounds_type(),
    const error_evaluator_type& error_evaluator = error_evaluator_type())
  {
    using value_type    = vector<scalar_type, 8>;
    using method_type   = explicit_method<tableau_type>;
    using problem_type  = initial_value_problem<scalar_type, mapped<value_type>, device_function<value_type(scalar_type, const value_type&)>>;
    using iterator_type = adaptive_step_iterator<method_type, problem_type, error_evaluator_type>;

    iterator_type iterator 
    {
      {
        lambda,                                                         // t0
        mapped<value_type>(ray.position.data()),                        // y0 - Valid as long as ray.position and ray.direction are contiguous in memory. 
        [&metric] __device__ (const scalar_type t, const value_type& y) // dy/dt = f(t,y)
        {
          value_type dydt;
          dydt.head(4) = y.tail(4);
          auto christoffel_symbols = metric.christoffel_symbols(y.head(4));
          for (auto i = 0; i < 4; ++i)
            for (auto j = 0; j < 4; ++j)
              for (auto k = 0; k < 4; ++k)
                dydt.tail(4)[k] -= christoffel_symbols(i, j, k) * y.tail(4)[i] * y.tail(4)[j];
          return dydt;
        }
      }, 
      lambda_step_size, 
      error_evaluator
    };
    
    for (std::size_t iteration = 0; iteration < iterations; ++iteration)
    {
      ++iterator;

      auto termination = metric.check_termination(ray.position, ray.direction);
      if (termination != termination_reason::none)
        return termination;
      if (!bounds.isEmpty() && !bounds.contains(ray.position))
        return termination_reason::out_of_bounds;
      if (ray.position.hasNaN() || ray.direction.hasNaN())
        return termination_reason::numeric_error;
    }
        
    return termination_reason::none;
  }
};
}
