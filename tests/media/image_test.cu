#include <doctest/doctest.h>

#include <astray/api.hpp>

TEST_CASE("ast::image")
{
  using pixel_type = ast::vector3<std::uint8_t>;
  using image_type = ast::image<pixel_type>;

  image_type image({256, 256}, pixel_type(0, 0, 0));
  for (auto x = 0; x < image.size[0]; ++x)
    for (auto y = 0; y < image.size[1]; ++y)
      image.at({x, y})          = pixel_type(static_cast<std::uint8_t>(x), static_cast<std::uint8_t>(y), 128);
  image.save("../data/outputs/tests/image_test.png");

  image_type other;
  other.load("../data/outputs/tests/image_test.png");
  for (auto x = 0; x < other.size[0]; ++x)
    for (auto y = 0; y < other.size[1]; ++y)
      REQUIRE(other.at({x, y}) == pixel_type(static_cast<std::uint8_t>(x), static_cast<std::uint8_t>(y), 128));
  
  REQUIRE(image.data == other.data);
  REQUIRE(image.size == other.size);
}