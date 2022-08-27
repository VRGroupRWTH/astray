#pragma once

#include <astray/parallel/thrust.hpp>

namespace ast
{
enum class shared_device_type
{
  cuda = THRUST_DEVICE_SYSTEM_CUDA,
  omp  = THRUST_DEVICE_SYSTEM_OMP ,
  tbb  = THRUST_DEVICE_SYSTEM_TBB ,
  cpp  = THRUST_DEVICE_SYSTEM_CPP
};

constexpr static shared_device_type shared_device = static_cast<shared_device_type>(THRUST_DEVICE_SYSTEM);
}