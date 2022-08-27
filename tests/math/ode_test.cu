#include <doctest/doctest.h>

#include <astray/api.hpp>

void test()
{
  using scalar_type             = float;
  using vector_type             = ast::vector3<scalar_type>;
  using tableau_type            = ast::dormand_prince_5_tableau<scalar_type>;
  using method_type             = ast::explicit_method<tableau_type>;
  using problem_type            = ast::initial_value_problem<scalar_type, vector_type>;
  using i_controller            = ast::integral_controller                        <scalar_type, tableau_type>;
  using pi_controller           = ast::proportional_integral_controller           <scalar_type, tableau_type>;
  using pid_controller          = ast::proportional_integral_derivative_controller<scalar_type, tableau_type>;
  using i_controller_iterator   = ast::adaptive_step_iterator<method_type, problem_type, i_controller  >;
  using pi_controller_iterator  = ast::adaptive_step_iterator<method_type, problem_type, pi_controller >;
  using pid_controller_iterator = ast::adaptive_step_iterator<method_type, problem_type, pid_controller>;
  
  std::vector<vector_type> input (1);
  std::vector<vector_type> output(1);

  thrust::device_vector<vector_type> device_data = input;
  thrust::for_each(device_data.begin(), device_data.end(), [ ] __device__ (auto& value)
  {
    constexpr auto sigma   = 10.0f;
    constexpr auto rho     = 28.0f;
    constexpr auto beta    = 8.0f / 3.0f;
    const     auto problem = problem_type
    {
      0.0f,                                     /* t0 */
      vector_type(16.0f, 16.0f, 16.0f),         /* y0 */
      [&] __device__ (const float t, const vector_type& y) /* y' = f(t, y) */
      {
        return vector_type(sigma * (y[1] - y[0]), y[0] * (rho - y[2]) - y[1], y[0] * y[1] - beta * y[2]); /* Lorenz system */
      }
    };
    
    auto iterator  = i_controller_iterator   {problem, 1.0f, {}};
    for (auto i = 0; i < 10000; ++i)
      ++iterator;
      
    auto iterator2 = pi_controller_iterator  {problem, 1.0f, {}};
    for (auto i = 0; i < 10000; ++i)
      ++iterator2;
      
    // auto iterator3 = pid_controller_iterator {problem, 1.0f, pid_controller()};
    // for (auto i = 0; i < 10000; ++i)
    //   ++iterator3;

    value = iterator.problem.value;
  });
  thrust::copy(device_data.begin(), device_data.end(), output.begin());
}

TEST_CASE("ast::ode")
{
  test();
}