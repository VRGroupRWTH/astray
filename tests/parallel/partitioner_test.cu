#include <doctest/doctest.h>

#include <astray/api.hpp>

TEST_CASE("ast::partitioner")
{
  using pixel_type       = ast::vector3<std::uint8_t>;
  using partitioner_type = ast::partitioner<2, std::int32_t, ast::vector2<std::int32_t>, true>;
  
  const ast::vector2<std::int32_t> domain_size = {640, 480};
  const std::int32_t               ranks       = 32;
  
  ast::image domain(domain_size, pixel_type(0, 0, 0));
  for (std::int32_t i = 0; i < ranks; ++i)
  {
    partitioner_type partitioner(i, ranks, domain_size);
    for (std::int32_t x = 0; x < partitioner.block_size()[0]; ++x)
      for (std::int32_t y = 0; y < partitioner.block_size()[1]; ++y)
        domain.at({partitioner.rank_offset()[0] + x, partitioner.rank_offset()[1] + y}) = pixel_type::Constant(static_cast<std::uint8_t>(255.0 * i / ranks));
  }
  domain.save("../data/outputs/tests/partitioner_test.jpg");
}