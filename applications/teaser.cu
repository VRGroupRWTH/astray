#include <cstdint>

#include <astray/api.hpp>

#include "common.hpp"

using scalar_type = double;

std::int32_t main(std::int32_t argc, char** argv)
{
  ast::mpi::environment environment;

  {
    settings_type<scalar_type, ast::metrics::minkowski<scalar_type>> settings;
    settings.image_size       = {3840, 2160};
    settings.iterations       = 40000;
    settings.lambda_step_size = static_cast<scalar_type>(0.001);
    settings.position         = {8, 0, 0};
    settings.background_image.load("../data/backgrounds/checkerboard_gray.png");
    std::get<ast::perspective_projection<scalar_type>>(settings.projection).fov_y = 100;

    const auto ray_tracer = make_ray_tracer(settings);
    const auto image      = ray_tracer->render_frame();
    if (ray_tracer->get_communicator().rank() == 0)
      image.save("../data/outputs/applications/teaser_minkowski.png");
  }
  {
    settings_type<scalar_type, ast::metrics::schwarzschild<scalar_type>> settings;
    settings.image_size       = {3840, 2160};
    settings.iterations       = 40000;
    settings.lambda_step_size = static_cast<scalar_type>(0.001);
    settings.position         = {8, 0, 0};
    settings.background_image.load("../data/backgrounds/checkerboard_gray.png");
    std::get<ast::perspective_projection<scalar_type>>(settings.projection).fov_y = 100;

    const auto ray_tracer = make_ray_tracer(settings);
    const auto image      = ray_tracer->render_frame();
    if (ray_tracer->get_communicator().rank() == 0)
      image.save("../data/outputs/applications/teaser_schwarzschild.png");
  }
  {
    settings_type<scalar_type, ast::metrics::kerr<scalar_type>> settings;
    settings.image_size       = {3840, 2160};
    settings.iterations       = 40000;
    settings.lambda_step_size = static_cast<scalar_type>(0.001);
    settings.position         = {8, 0, 0};
    settings.background_image.load("../data/backgrounds/checkerboard_gray.png");
    std::get<ast::perspective_projection<scalar_type>>(settings.projection).fov_y = 100;

    const auto ray_tracer = make_ray_tracer(settings);
    const auto image      = ray_tracer->render_frame();
    if (ray_tracer->get_communicator().rank() == 0)
      image.save("../data/outputs/applications/teaser_kerr.png");
  }
  {
    settings_type<scalar_type, ast::metrics::reissner_nordstroem<scalar_type>> settings;
    settings.image_size       = {3840, 2160};
    settings.iterations       = 40000;
    settings.lambda_step_size = static_cast<scalar_type>(0.001);
    settings.position         = {8, 0, 0};
    settings.background_image.load("../data/backgrounds/checkerboard_gray.png");
    std::get<ast::perspective_projection<scalar_type>>(settings.projection).fov_y = 100;

    const auto ray_tracer = make_ray_tracer(settings);
    const auto image      = ray_tracer->render_frame();
    if (ray_tracer->get_communicator().rank() == 0)
      image.save("../data/outputs/applications/teaser_reissner_nordstroem.png");
  }
  {
    settings_type<scalar_type, ast::metrics::morris_thorne<scalar_type>> settings;
    settings.image_size       = {3840, 2160};
    settings.iterations       = 40000;
    settings.lambda_step_size = static_cast<scalar_type>(0.001);
    settings.position         = {1, 0, 0};
    settings.background_image.load("../data/backgrounds/checkerboard_gray.png");
    std::get<ast::perspective_projection<scalar_type>>(settings.projection).fov_y = 100;

    const auto ray_tracer = make_ray_tracer(settings);
    const auto image      = ray_tracer->render_frame();
    if (ray_tracer->get_communicator().rank() == 0)
      image.save("../data/outputs/applications/teaser_morris_thorne.png");
  }
  {
    settings_type<scalar_type, ast::metrics::kastor_traschen<scalar_type>> settings;
    settings.image_size       = {3840, 2160};
    settings.iterations       = 40000;
    settings.lambda_step_size = static_cast<scalar_type>(0.001);
    settings.position         = {5, 0, 0};
    settings.background_image.load("../data/backgrounds/checkerboard_gray.png");
    std::get<ast::perspective_projection<scalar_type>>(settings.projection).fov_y = 100;

    const auto ray_tracer = make_ray_tracer(settings);
    const auto image      = ray_tracer->render_frame();
    if (ray_tracer->get_communicator().rank() == 0)
      image.save("../data/outputs/applications/teaser_kastor_traschen.png");
  }

  return 0;
}