#include <doctest/doctest.h>

#include <astray/api.hpp>

TEST_CASE("ast::mpi")
{
  if constexpr (ast::distributed_device == ast::distributed_device_type::mpi)
  {
    ast::mpi::environment  environment ;
    ast::mpi::communicator communicator;
    
    std::int32_t value = 0;
    if (communicator.rank() == 0)
      value = 42;

    MPI_Bcast(&value, 1, MPI_INT, 0, communicator.native());
    REQUIRE(value == 42);
  }
}