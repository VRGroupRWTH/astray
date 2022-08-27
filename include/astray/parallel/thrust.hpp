#pragma once

#include <astray/utility/warning.hpp>

DISABLE_WARNING_PUSH
DISABLE_WARNING_PRAGMAS
DISABLE_WARNING_NO_DEPRECATED_GPU_TARGETS
#include <thrust/iterator/counting_iterator.h>
#include <thrust/device_vector.h>
#include <thrust/for_each.h>

#ifdef __CUDACC__
#include <nvfunctional>
namespace ast
{
template <typename type>
using device_function = nvstd::function<type>;
}
#else
#include <functional>
namespace ast
{
template <typename type>
using device_function = std::function<type>;
}
#endif

DISABLE_WARNING_POP

#ifndef __constant__
#define __constant__
#endif
