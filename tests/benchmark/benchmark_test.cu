#include <doctest/doctest.h>

#include <numeric>
#include <vector>

#include <astray/api.hpp>

TEST_CASE("ast::benchmark")
{
  auto record = ast::benchmark<double, std::milli, std::chrono::high_resolution_clock> ([ ] ()
  {
    std::vector<double> vector(static_cast<std::size_t>(1e8)); // 0.8 GB.
    std::iota(vector.begin(), vector.end(), 0);
    vector.clear();
  }, static_cast<std::size_t>(10));
  record.to_csv("../data/outputs/tests/benchmark_test_micro.csv");

  auto session = ast::benchmark        <double, std::milli, std::chrono::high_resolution_clock> ([ ] (auto& recorder)
  {
    std::vector<double> vector;
    recorder.record("resize"  , [&]
    {
      vector.resize(static_cast<std::size_t>(1e8)); // 0.8 GB.
    });
    recorder.record("generate", [&]
    {
      std::iota(vector.begin(), vector.end(), 0);
    });
    recorder.record("clear"   , [&]
    {
      vector.clear();
    });
  }, static_cast<std::size_t>(10));
  session.to_csv("../data/outputs/tests/benchmark_test_macro.csv");

  if constexpr (ast::distributed_device == ast::distributed_device_type::mpi)
  {
    ast::mpi::environment  environment ;
    ast::mpi::communicator communicator;

    auto mpi_session = ast::benchmark_mpi<double, std::milli, std::chrono::high_resolution_clock> ([ ] (auto& recorder)
    {
      std::vector<double> vector;
      recorder.record("resize"  , [&]
      {
        vector.resize(static_cast<std::size_t>(1e8)); // 0.8 GB.
      });
      recorder.record("generate", [&]
      {
        std::iota(vector.begin(), vector.end(), 0);
      });
      recorder.record("clear"   , [&]
      {
        vector.clear();
      });
    }, static_cast<std::size_t>(10));

    if (communicator.rank() == 0)
      mpi_session.to_csv("../data/outputs/tests/benchmark_test_mpi.csv");
  }
}