#include <doctest/doctest.h>

#include <astray/api.hpp>

TEST_CASE("ast::ray_tracer")
{
  using scalar_type     = float;
  using metric_type     = ast::metrics::kerr<scalar_type>;
  using geodesic_type   = ast::geodesic<scalar_type, ast::runge_kutta_4_tableau<scalar_type>>;
  using ray_tracer_type = ast::ray_tracer<metric_type, geodesic_type>;

  auto ray_tracer = ray_tracer_type
  (
    {640, 480},
    {},
    1000,
    0.01f,
    0.0f,
    {},
    {},
    false
  );

  ray_tracer.set_background(ray_tracer_type::image_type("../data/backgrounds/checkerboard.png"));
  ray_tracer.get_observer().get_transform().translation = {0.1f, 0.1f, 10.0f};
  ray_tracer.get_observer().get_transform().look_at({0.0f, 0.0f, 0.0f});
  ray_tracer.get_observer().set_projection(ast::perspective_projection<scalar_type> {ast::to_radians(75.0f), 640.0f / 480.0f});

  const auto image = ray_tracer.render_frame();
  if (ray_tracer.get_communicator().rank() == 0)
    image.save("../data/outputs/tests/ray_tracer_test.jpg");
  for (auto i = 1; i < ray_tracer.get_communicator().size(); ++i)
    if (ray_tracer.get_communicator().rank() == i)
      image.save("../data/outputs/tests/ray_tracer_test_rank" + std::to_string(i) + ".jpg");
}