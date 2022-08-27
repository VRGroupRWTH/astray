#pragma once

#define _USE_MATH_DEFINES

#include <cmath>
#include <math.h>

namespace ast
{
template <typename type>
constexpr type to_radians_multiplier = static_cast<type>(M_PI / 180.0);
template <typename type>
constexpr type to_degrees_multiplier = static_cast<type>(180.0 / M_PI);

template <typename type>
constexpr type to_radians(const type degrees)
{
  return degrees * to_radians_multiplier<type>;
}
template <typename type>
constexpr type to_degrees(const type radians)
{
  return radians * to_degrees_multiplier<type>;
}

template <typename type>
constexpr void wrap_angle(type& value, const type maximum, const bool closed_interval = false)
{
  if (value == 0)
    return;
  
  value = std::fmod(value, maximum);

  if (closed_interval ? value <= 0 : value < 0)
    value += maximum;
}
}