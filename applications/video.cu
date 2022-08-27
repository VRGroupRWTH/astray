#include <cstdint>
#include <iostream>

#include <astray/api.hpp>

#include "common.hpp"

using scalar_type = float;

std::int32_t main(std::int32_t argc, char** argv)
{
  settings_type<scalar_type, ast::metrics::schwarzschild<scalar_type>> settings;
  settings.image_size = {960, 540};
  settings.position   = {5, 0, -5};
  settings.background_image.load("../data/backgrounds/checkerboard_gray.png");

  const auto ray_tracer = make_ray_tracer(settings);

  std::optional<ast::video> video(std::nullopt);
  if (ray_tracer->get_communicator().rank() == 0)
    video.emplace("../data/outputs/applications/video.mp4", ray_tracer->get_image_size(), 60);
  
  constexpr auto frames(1000);
  for (auto i = 0; i < frames; ++i)
  {
    if (i % 10 == 0)
      std::cout << i << "/" << frames - 1 << "\n";

    auto image = ray_tracer->render_frame();
    if (ray_tracer->get_communicator().rank() == 0)
      video->append(image);

    ray_tracer->get_observer().get_transform().translation[2] += 0.01f;
    ray_tracer->get_observer().get_transform().look_at({0, 0, 0});
  }

  return 0;
}