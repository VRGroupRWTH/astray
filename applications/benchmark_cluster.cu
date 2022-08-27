#include <cstdint>
#include <iostream>

#include <astray/api.hpp>

#include "common.hpp"

using scalar_type = float;

template <typename scalar_type, typename metric_type, typename motion_type>
constexpr auto run_benchmark  (
  settings_type<scalar_type, metric_type, motion_type> settings   , 
  std::size_t                                          repeats    , 
  const std::string&                                   device_name,
  const std::string&                                   metric_name)
{
  using ray_tracer_type = ast::ray_tracer<metric_type, motion_type>;
  using image_type      = typename ray_tracer_type::image_type;
  using image_size_type = typename ray_tracer_type::image_size_type;

  std::get<ast::perspective_projection<scalar_type>>(settings.projection).aspect_ratio = static_cast<scalar_type>(1);

  auto              ray_tracer = make_ray_tracer(settings);
  const image_type* image;
  
  auto iteration  = 0;
  auto image_size = image_size_type();
#ifdef ASTRAY_USE_MPI
  auto session    = ast::benchmark_mpi<scalar_type, std::milli, std::chrono::high_resolution_clock>([&] (auto& recorder)
#else
  auto session    = ast::benchmark    <scalar_type, std::milli, std::chrono::high_resolution_clock>([&] (auto& recorder)
#endif
  {
    image_size = image_size_type(512, 512);
    ray_tracer->set_image_size(image_size);
    if (ray_tracer->get_communicator().rank() == 0)
      std::cout << "Recording. Processes: " << ray_tracer->get_communicator().size() << ". Iteration: " << iteration << ". Device: " << device_name << ". Metric: " << metric_name << ". Image Size: " << image_size[0] << "x" << image_size[1] << ".\n";
    recorder.record(metric_name + "," + std::to_string(image_size[0]) + "," + std::to_string(image_size[1]), [&]
    {
      image = &ray_tracer->render_frame();
    });

    image_size = image_size_type(724, 724);
    ray_tracer->set_image_size(image_size);
    if (ray_tracer->get_communicator().rank() == 0)
      std::cout << "Recording. Processes: " << ray_tracer->get_communicator().size() << ". Iteration: " << iteration << ". Device: " << device_name << ". Metric: " << metric_name << ". Image Size: " << image_size[0] << "x" << image_size[1] << ".\n";
    recorder.record(metric_name + "," + std::to_string(image_size[0]) + "," + std::to_string(image_size[1]), [&]
    {
      image = &ray_tracer->render_frame();
    });

    image_size = image_size_type(1024, 1024);
    ray_tracer->set_image_size(image_size);
    if (ray_tracer->get_communicator().rank() == 0)
      std::cout << "Recording. Processes: " << ray_tracer->get_communicator().size() << ". Iteration: " << iteration << ". Device: " << device_name << ". Metric: " << metric_name << ". Image Size: " << image_size[0] << "x" << image_size[1] << ".\n";
    recorder.record(metric_name + "," + std::to_string(image_size[0]) + "," + std::to_string(image_size[1]), [&]
    {
      image = &ray_tracer->render_frame();
    });

    image_size = image_size_type(1448, 1448);
    ray_tracer->set_image_size(image_size);
    if (ray_tracer->get_communicator().rank() == 0)
      std::cout << "Recording. Processes: " << ray_tracer->get_communicator().size() << ". Iteration: " << iteration << ". Device: " << device_name << ". Metric: " << metric_name << ". Image Size: " << image_size[0] << "x" << image_size[1] << ".\n";
    recorder.record(metric_name + "," + std::to_string(image_size[0]) + "," + std::to_string(image_size[1]), [&]
    {
      image = &ray_tracer->render_frame();
    });

    image_size = image_size_type(2048, 2048);
    ray_tracer->set_image_size(image_size);
    if (ray_tracer->get_communicator().rank() == 0)
      std::cout << "Recording. Processes: " << ray_tracer->get_communicator().size() << ". Iteration: " << iteration << ". Device: " << device_name << ". Metric: " << metric_name << ". Image Size: " << image_size[0] << "x" << image_size[1] << ".\n";
    recorder.record(metric_name + "," + std::to_string(image_size[0]) + "," + std::to_string(image_size[1]), [&]
    {
      image = &ray_tracer->render_frame();
    });

    //image_size = image_size_type(2896, 2896);
    //ray_tracer->set_image_size(image_size);
    //if (ray_tracer->get_communicator().rank() == 0)
    //  std::cout << "Recording. Processes: " << ray_tracer->get_communicator().size() << ". Iteration: " << iteration << ". Device: " << device_name << ". Metric: " << metric_name << ". Image Size: " << image_size[0] << "x" << image_size[1] << ".\n";
    //recorder.record(metric_name + "," + std::to_string(image_size[0]) + "," + std::to_string(image_size[1]), [&]
    //{
    //  image = &ray_tracer->render_frame();
    //});
    //
    //image_size = image_size_type(4096, 4096);
    //if (ray_tracer->get_communicator().rank() == 0)
    //  std::cout << "Recording. Processes: " << ray_tracer->get_communicator().size() << ". Iteration: " << iteration << ". Device: " << device_name << ". Metric: " << metric_name << ". Image Size: " << image_size[0] << "x" << image_size[1] << ".\n";
    //recorder.record(metric_name + "," + std::to_string(image_size[0]) + "," + std::to_string(image_size[1]), [&]
    //{
    //  image = &ray_tracer->render_frame();
    //});

    ++iteration;
  }, repeats);

  if (ray_tracer->get_communicator().rank() == 0)
    image->save("./benchmark_cluster_" + device_name + "_" + metric_name + ".png"); // Just put it in the execution directory.

  return session;
}

std::int32_t main(std::int32_t argc, char** argv)
{
  ast::mpi::environment  environment ;
  ast::mpi::communicator communicator;
  
  constexpr auto runs = 5;
  
  std::string device_name;
  if      constexpr (ast::shared_device == ast::shared_device_type::cpp )
    device_name = "cpp" ;
  else if constexpr (ast::shared_device == ast::shared_device_type::cuda)
    device_name = "cuda";
  else if constexpr (ast::shared_device == ast::shared_device_type::omp )
    device_name = "omp" ;
  else if constexpr (ast::shared_device == ast::shared_device_type::tbb )
    device_name = "tbb" ;
  
  auto benchmark1 = run_benchmark(settings_type<scalar_type, ast::metrics::minkowski          <scalar_type>>(), runs, device_name, "minkowski"          );
  auto benchmark2 = run_benchmark(settings_type<scalar_type, ast::metrics::schwarzschild      <scalar_type>>(), runs, device_name, "schwarzschild"      );
  auto benchmark3 = run_benchmark(settings_type<scalar_type, ast::metrics::kerr               <scalar_type>>(), runs, device_name, "kerr"               );
  //auto benchmark4 = run_benchmark(settings_type<scalar_type, ast::metrics::reissner_nordstroem<scalar_type>>(), runs, device_name, "reissner_nordstroem");
  //auto benchmark5 = run_benchmark(settings_type<scalar_type, ast::metrics::morris_thorne      <scalar_type>>(), runs, device_name, "morris_thorne"      );
  auto benchmark6 = run_benchmark(settings_type<scalar_type, ast::metrics::kastor_traschen    <scalar_type>>(), runs, device_name, "kastor_traschen"    );

  if (communicator.rank() == 0)
  {
    std::ofstream stream("../data/outputs/performance/benchmark_cluster_" + device_name + "_" + std::to_string(communicator.size()) + ".csv");
    stream << "rank,metric,width,height,";
    for (auto i = 0; i < runs; ++i)
      stream << "run_" << i << ",";
    stream << "mean,variance,standard deviation\n";
    stream << benchmark1.to_string();
    stream << benchmark2.to_string();
    stream << benchmark3.to_string();
    //stream << benchmark4.to_string();
    //stream << benchmark5.to_string();
    stream << benchmark6.to_string();
  }

  return 0;
}
