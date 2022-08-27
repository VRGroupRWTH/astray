#include <astray/api.hpp>

int main(int argc, char** argv)
{
  using scalar_type   = double;
  using metric_type   = ast::metrics::schwarzschild<scalar_type>;
  using geodesic_type = ast::geodesic<scalar_type, ast::runge_kutta_4_tableau<scalar_type>>;

  ast::ray_tracer<metric_type, geodesic_type> ray_tracer;
  ray_tracer.set_image_size      ({1920, 1080});
  ray_tracer.set_lambda_step_size(0.01);
  ray_tracer.set_iterations      (1000);

  auto& transform       = ray_tracer.get_observer().get_transform();
  transform.translation = {0.0, 0.0, -10.0};
  transform.look_at      ({0.0, 0.0,   0.0});

  const auto image = ray_tracer.render_frame();
  image.save("example.png");
}