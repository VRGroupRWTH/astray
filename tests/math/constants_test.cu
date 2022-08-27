#include <doctest/doctest.h>

#include <astray/api.hpp>

void add_pi(const std::vector<double>& input, std::vector<double>& output)
{
  thrust::device_vector<double> device_data = input;
  thrust::for_each(device_data.begin(), device_data.end(), [ ] __device__ (auto& value)
  {
    value += ast::constants<double>::pi;
  });
  thrust::copy(device_data.begin(), device_data.end(), output.begin());
}

TEST_CASE("ast::constants")
{
  std::vector input      {1.0};
  std::vector output     {0.0};
  std::vector comparison {1.0 + M_PI};

  add_pi(input, output);

  REQUIRE(output == comparison);
}