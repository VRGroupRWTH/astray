#include <doctest/doctest.h>

#include <astray/api.hpp>

TEST_CASE("ast::observer")
{
  using scalar_type     = double;
  using observer_type   = ast::observer<scalar_type>;
  using ray_type        = observer_type::ray_type;
  using pixel_type      = ast::vector3<std::uint8_t>;
  using image_type      = ast::image<pixel_type>;
  using image_size_type = image_type::size_type;

  image_size_type size {640, 480};

  {
    observer_type observer;
    observer.set_projection(ast::perspective_projection<scalar_type>{ast::to_radians(75.0), static_cast<scalar_type>(size[0]) / size[1]});

    auto& rays      = observer.generate_rays(size, size, image_size_type::Zero().eval());
    auto  rays_host = std::vector<ray_type>(size.prod());
    thrust::copy(rays.begin(), rays.end(), rays_host.begin());

    auto image      = image_type({size});

    for (auto x = 0; x < size[0]; ++x)
      for (auto y = 0; y < size[1]; ++y)
      {
        auto& ray = rays_host[ast::ravel_multi_index<image_size_type, true>({x, y}, size)];
        image.at({x, y}) = (ray.position.tail(3).cwiseAbs() * 255).cast<std::uint8_t>();
      }
    image.save("../data/outputs/tests/observer_test_perspective_positions.jpg");
    
    for (auto x = 0; x < size[0]; ++x)
      for (auto y = 0; y < size[1]; ++y)
      {
        auto& ray = rays_host[ast::ravel_multi_index<image_size_type, true>({x, y}, size)];
        image.at({x, y}) = (ray.direction.tail(3).cwiseAbs() * 255).cast<std::uint8_t>();
      }
    image.save("../data/outputs/tests/observer_test_perspective_directions.jpg");
  }

  {
    observer_type observer;
    observer.set_projection(ast::orthographic_projection<scalar_type>{static_cast<scalar_type>(1), static_cast<scalar_type>(size[0]) / size[1]});

    auto& rays      = observer.generate_rays(size, size, image_size_type::Zero().eval());
    auto  rays_host = std::vector<ray_type>(size.prod());
    thrust::copy(rays.begin(), rays.end(), rays_host.begin());

    auto image      = image_type({size});

    for (auto x = 0; x < size[0]; ++x)
      for (auto y = 0; y < size[1]; ++y)
      {
        auto& ray = rays_host[ast::ravel_multi_index<image_size_type, true>({x, y}, size)];
        image.at({x, y}) = (ray.position.tail(3).cwiseAbs() * 255).cast<std::uint8_t>();
      }
    image.save("../data/outputs/tests/observer_test_orthographic_positions.jpg");
    
    for (auto x = 0; x < size[0]; ++x)
      for (auto y = 0; y < size[1]; ++y)
      {
        auto& ray = rays_host[ast::ravel_multi_index<image_size_type, true>({x, y}, size)];
        image.at({x, y}) = (ray.direction.tail(3).cwiseAbs() * 255).cast<std::uint8_t>();
      }
    image.save("../data/outputs/tests/observer_test_orthographic_directions.jpg");
  }
}