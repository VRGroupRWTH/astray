![Teaser](./data/docs/teaser.png)

### Astray
A performance-portable geodesic ray tracing library.

### Requirements
- CMake 3.23+
- Optional: CUDA 11.6+
- Optional: MPI 2.0+

### Installation
- Clone the repository.
- Run `bootstrap.[bat|sh]`. This will install all dependencies except CUDA and MPI, and create the project under the `./build` directory.
- Optional: Run cmake on the `./build` directory.
  - Toggle `ASTRAY_BUILD_APPLICATIONS` for building the example applications.
  - Toggle `ASTRAY_BUILD_TESTS` for building the tests.
  - Toggle `ASTRAY_USE_FFMPEG` for video support.
  - Toggle `ASTRAY_USE_MPI` for MPI support. The build will ask for the location of MPI upon enabling this option.
  - Toggle `ASTRAY_DEVICE_SYSTEM` for CUDA/OMP/TBB support. The build will ask for the location of the respective library upon enabling this option.
  - Remember to generate or run `bootstrap.bat` after changes. You can ignore the cmake developer errors as long as generation is successful.

### Example Usage: Rendering a Schwarzschild Black Hole
```cpp
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
```
See the tests for more.

### UI Application
See [Astrid](https://github.com/VRGroupRWTH/astrid) for a UI application built on Astray.

### Future Work
- Fermi Walker transport.
- Precompute the Christoffel symbols into a 4D image.
- More tests.