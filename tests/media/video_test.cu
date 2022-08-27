#include <doctest/doctest.h>

#include <astray/api.hpp>

TEST_CASE("ast::video")
{
  using pixel_type = ast::vector3<std::uint8_t>;
  using frame_type = ast::image<pixel_type>;

  ast::video video("../data/outputs/tests/video_test.mp4", {256, 256}, 60);
  frame_type frame({256, 256}, pixel_type(0, 0, 0));

  for (auto i = 0; i < 256; ++i)
  {
    for (auto x = 0; x < frame.size[0]; ++x)
      for (auto y = 0; y < frame.size[1]; ++y)
        frame.at({x, y}) = pixel_type(static_cast<std::uint8_t>(x), static_cast<std::uint8_t>(y), static_cast<std::uint8_t>(i));
    video.append(frame);
  }
}