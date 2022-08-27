#include <doctest/doctest.h>

#include <astray/api.hpp>

template <typename type>
void add_one(const std::vector<type>& input, std::vector<type>& output)
{
  thrust::device_vector<type> device_data = input;
  thrust::for_each(device_data.begin(), device_data.end(), [ ] __device__ (auto& value)
  {
    value += type::Ones();
  });
  thrust::copy(device_data.begin(), device_data.end(), output.begin());
}

TEST_CASE("ast::linear_algebra")
{
  {
    using type = ast::vector3<double>;

    std::vector<type> input      { type::Ones    ()  };
    std::vector<type> output     { type::Zero    ()  };
    std::vector<type> comparison { type::Constant(2) };
    add_one(input, output);
    REQUIRE(output == comparison);
  }
  
  {
    using type = ast::matrix33<double>;

    std::vector<type> input      { type::Ones    ()  };
    std::vector<type> output     { type::Zero    ()  };
    std::vector<type> comparison { type::Constant(2) };
    add_one(input, output);
    REQUIRE(output == comparison);
  }

  // TODO: Tests for tensors, AABBs, quaternions, angle axes.
}