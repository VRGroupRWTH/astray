#pragma once

namespace ast
{
enum class distributed_device_type
{
  none,
  mpi
};

constexpr static distributed_device_type distributed_device = 
#ifdef ASTRAY_USE_MPI 
  distributed_device_type::mpi;
#else
  distributed_device_type::none;
#endif
}