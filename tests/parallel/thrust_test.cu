#include <doctest/doctest.h>

#include <astray/api.hpp>

void add_one(const std::vector<double>& input, std::vector<double>& output)
{
  thrust::device_vector<double> device_data = input;
  thrust::for_each(device_data.begin(), device_data.end(), [ ] __device__ (auto& value)
  {
    value += 1.0;
  });
  thrust::copy(device_data.begin(), device_data.end(), output.begin());
}

TEST_CASE("ast::thrust")
{
  std::vector<double> input     (1000);
  std::vector<double> output    (1000);
  std::vector<double> comparison(1000);

  std::iota(input     .begin(), input     .end(), 0.0);
  std::iota(comparison.begin(), comparison.end(), 1.0);

  add_one(input, output);

  REQUIRE(output == comparison);
}