#include <doctest/doctest.h>

#include <astray/api.hpp>

void test()
{
  using scalar_type           = double;
  using vector_type           = ast::vector4<scalar_type>;
  using ray_type              = ast::ray    <vector_type>;

  using tableau_type          = ast::dormand_prince_5_tableau<scalar_type>;
  using error_controller_type = ast::integral_controller<scalar_type, tableau_type>;
  using geodesic_type         = ast::geodesic<scalar_type, tableau_type, error_controller_type>;

  std::vector<ray_type> ray{{vector_type(0, 0, 0, 1), vector_type(0, 0, 0, 1)}};

  thrust::device_vector<ray_type> device_ray = ray;
  thrust::for_each(device_ray.begin(), device_ray.end(), [ ] __device__ (auto& value)
  {
    ast::metrics::minkowski<scalar_type> metric;
    geodesic_type::integrate(value, metric, 10, 0.1);
  });
  thrust::copy(device_ray.begin(), device_ray.end(), ray.begin());
}

TEST_CASE("ast::geodesic")
{
  test();

  // TODO
}