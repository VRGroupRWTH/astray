#pragma once

#include <cstddef>

#include <astray/api.hpp>

template <
  typename scalar_type = float,
  typename metric_type = ast::metrics::kerr<scalar_type>, 
  typename motion_type = ast::geodesic<scalar_type, ast::runge_kutta_4_tableau<scalar_type>>>
struct settings_type
{
  using ray_tracer_type      = ast::ray_tracer<metric_type, motion_type>;
  using observer_type        = typename ray_tracer_type::observer_type;
  using image_type           = typename ray_tracer_type::image_type;
  using image_size_type      = typename ray_tracer_type::image_size_type;
  using bounds_type          = typename ray_tracer_type::bounds_type;
  using error_evaluator_type = typename ray_tracer_type::error_evaluator_type;
  using transform_type       = typename observer_type  ::transform_type;
  using projection_type      = typename observer_type  ::projection_type;
  using vector_type          = typename transform_type ::vector_type;
  
  image_size_type      image_size       = {1920, 1080};
  metric_type          metric           = {};
  std::size_t          iterations       = 1000;
  scalar_type          lambda_step_size = static_cast<scalar_type>(0.01);
  scalar_type          lambda           = static_cast<scalar_type>(0);
  bounds_type          bounds           = {};
  error_evaluator_type error_evaluator  = {};
  bool                 debug            = false;

  vector_type          position         = vector_type(5, 0, 0);
  vector_type          rotation         = vector_type(0, 0, 0);
  bool                 look_at_origin   = true;
  scalar_type          coordinate_time  = static_cast<scalar_type>(0);
  projection_type      projection       = ast::perspective_projection<scalar_type> {ast::to_radians<scalar_type>(120), static_cast<scalar_type>(image_size[0]) / image_size[1]};
  image_type           background_image = image_type();
};

template <typename scalar_type, typename metric_type, typename motion_type>
constexpr auto make_ray_tracer(const settings_type<scalar_type, metric_type, motion_type>& settings)
{
#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
  cudaDeviceReset();
#endif
  
  auto ray_tracer = std::make_unique<ast::ray_tracer<metric_type, motion_type>>(
    settings.image_size      ,
    settings.metric          ,
    settings.iterations      ,
    settings.lambda_step_size,
    settings.lambda          ,
    settings.bounds          ,
    settings.error_evaluator ,
    settings.debug           );
  ray_tracer->get_observer().get_transform().translation = settings.position;
  ray_tracer->get_observer().get_transform().rotation_from_euler(settings.rotation);
  if (settings.look_at_origin)
    ray_tracer->get_observer().get_transform().look_at(ast::transform<scalar_type>::vector_type::Zero());
  ray_tracer->get_observer().set_coordinate_time(settings.coordinate_time);
  ray_tracer->get_observer().set_projection     (settings.projection);
  ray_tracer->set_background(settings.background_image);
  return ray_tracer;
}