#pragma once

#include <cmath>

#include <astray/math/angle.hpp>
#include <astray/math/constants.hpp>
#include <astray/math/linear_algebra.hpp>

namespace ast
{
enum class coordinate_system_type
{
  boyer_lindquist   , // t [-INF, INF], r     [   0, INF], theta [   0,  PI], phi [0   , 2PI)
  cartesian         , // t [-INF, INF], x     [-INF, INF], y     [-INF, INF], z   [-INF, INF]
  cylindrical       , // t [-INF, INF], rho   [   0, INF], phi   [   0, 2PI), z   [-INF, INF]
  prolate_spheroidal, // t [-INF, INF], sigma [   1, INF], tau   [  -1,   1], phi [0   , 2PI)
  spherical           // t [-INF, INF], r     [   0, INF], theta [   0,  PI], phi [0   , 2PI)
};

template <coordinate_system_type system, typename type>
__device__ __host__ constexpr void wrap_angles(type& value)
{
  if      constexpr (system == coordinate_system_type::boyer_lindquist)
  {
    wrap_angle(value[2], constants<typename type::value_type>::pi    , true );
    wrap_angle(value[3], constants<typename type::value_type>::two_pi, false);
  }
  else if constexpr (system == coordinate_system_type::cylindrical)
  {
    wrap_angle(value[2], constants<typename type::value_type>::two_pi, false);
  }
  else if constexpr (system == coordinate_system_type::prolate_spheroidal)
  {
    wrap_angle(value[3], constants<typename type::value_type>::two_pi, false);
  }
  else if constexpr (system == coordinate_system_type::spherical)
  {
    wrap_angle(value[2], constants<typename type::value_type>::pi    , true );
    wrap_angle(value[3], constants<typename type::value_type>::two_pi, false);
  }
}

template<coordinate_system_type source, coordinate_system_type target, typename type, typename scalar = typename type::value_type>
__device__ __host__ constexpr void convert    (type& value)
{
  if        constexpr (source == coordinate_system_type::cartesian)
  {
    if      constexpr (target == coordinate_system_type::cylindrical)
    {
      const scalar rho = std::sqrt (static_cast<scalar>(std::pow(value[1], 2)) + static_cast<scalar>(std::pow(value[2], 2)));
      const scalar phi = std::atan2(value[2], value[1]);

      value[1] = rho;
      value[2] = phi;

      wrap_angles<target>(value);
    }
    else if constexpr (target == coordinate_system_type::spherical)
    {
      const scalar r     = std::sqrt (
        static_cast<scalar>(std::pow(value[1], 2)) + 
        static_cast<scalar>(std::pow(value[2], 2)) + 
        static_cast<scalar>(std::pow(value[3], 2)));
      const scalar theta = std::acos (value[3] / r);
      const scalar phi   = std::atan2(value[2], value[1]);

      value[1] = r;
      value[2] = theta;
      value[3] = phi;

      wrap_angles<target>(value);
    }
  }
  else if   constexpr (source == coordinate_system_type::cylindrical)
  {
    if      constexpr (target == coordinate_system_type::cartesian)
    {
      wrap_angles<source>(value);

      const scalar x = value[1] * std::cos(value[2]);
      const scalar y = value[1] * std::sin(value[2]);

      value[1] = x;
      value[2] = y;
    }
    else if constexpr (target == coordinate_system_type::spherical)
    {
      wrap_angles<source>(value);

      const scalar r     = std::sqrt(
        static_cast<scalar>(std::pow(value[1], 2)) +
        static_cast<scalar>(std::pow(value[3], 2)));
      const scalar theta = std::atan2(value[1], value[3]);
      const scalar phi   = value[2];

      value[1] = r;
      value[2] = theta;
      value[3] = phi;

      wrap_angles<target>(value);
    }
  }
  else if   constexpr (source == coordinate_system_type::spherical)
  {
    if      constexpr (target == coordinate_system_type::cartesian)
    {
      wrap_angles<source>(value);

      const scalar x = value[1] * std::sin(value[2]) * std::cos(value[3]);
      const scalar y = value[1] * std::sin(value[2]) * std::sin(value[3]);
      const scalar z = value[1] * std::cos(value[2]);

      value[1] = x;
      value[2] = y;
      value[3] = z;
    }
    else if constexpr (target == coordinate_system_type::cylindrical)
    {
      wrap_angles<source>(value);

      const scalar rho = value[1] * std::sin(value[2]);
      const scalar phi = value[3];
      const scalar z   = value[1] * std::cos(value[2]);

      value[1] = rho;
      value[2] = phi;
      value[3] = z;

      wrap_angles<target>(value);
    }
  }
}
template<coordinate_system_type source, coordinate_system_type target, typename type, typename scalar = typename type::value_type>
__device__ __host__ constexpr void convert    (type& value, const scalar free_parameter)
{
  if        constexpr (source == coordinate_system_type::boyer_lindquist)
  {
    if      constexpr (target == coordinate_system_type::cartesian)
    {
      wrap_angles<source>(value);

      const scalar sq_r2_p_a2 = std::sqrt(
        static_cast<scalar>(std::pow(value[1]      , 2)) + 
        static_cast<scalar>(std::pow(free_parameter, 2)));

      const scalar x = sq_r2_p_a2 * std::sin(value[2]) * std::cos(value[3]);
      const scalar y = sq_r2_p_a2 * std::sin(value[2]) * std::sin(value[3]);
      const scalar z = value[1] * std::cos(value[2]);

      value[1] = x;
      value[2] = y;
      value[3] = z;
    }
    else if constexpr (target == coordinate_system_type::cylindrical)
    {
      convert<source, coordinate_system_type::cartesian>(value, free_parameter);
      convert<coordinate_system_type::cartesian, target>(value);
    }
    else if constexpr (target == coordinate_system_type::spherical)
    {
      convert<source, coordinate_system_type::cartesian>(value, free_parameter);
      convert<coordinate_system_type::cartesian, target>(value);
    }
  }
  else if   constexpr (source == coordinate_system_type::cartesian)
  {
    if      constexpr (target == coordinate_system_type::boyer_lindquist)
    {
      const scalar w = 
        static_cast<scalar>(std::pow(value[1]      , 2)) + 
        static_cast<scalar>(std::pow(value[2]      , 2)) + 
        static_cast<scalar>(std::pow(value[3]      , 2)) - 
        static_cast<scalar>(std::pow(free_parameter, 2));

      const scalar r     = std::sqrt (static_cast<scalar>(0.5) * (w + std::sqrt(static_cast<scalar>(std::pow(w, 2)) + static_cast<scalar>(4) * static_cast<scalar>(std::pow(free_parameter, 2)) * static_cast<scalar>(std::pow(value[3], 2)))));
      const scalar theta = std::acos (value[3] / r);
      const scalar phi   = std::atan2(value[2], value[1]);

      value[1] = r;
      value[2] = theta;
      value[3] = phi;

      wrap_angles<target>(value);
    }
    else if constexpr (target == coordinate_system_type::prolate_spheroidal)
    {
      const scalar x_squared        = static_cast<scalar>(std::pow(value[1], 2));
      const scalar y_squared        = static_cast<scalar>(std::pow(value[2], 2));
      const scalar first_component  = std::sqrt(x_squared + y_squared + static_cast<scalar>(std::pow(value[3] + free_parameter, 2)));
      const scalar second_component = std::sqrt(x_squared + y_squared + static_cast<scalar>(std::pow(value[3] - free_parameter, 2)));
      const scalar _2a              = 2 * free_parameter;

      const scalar sigma = (first_component + second_component) / _2a;
      const scalar tau   = (first_component - second_component) / _2a; 
      const scalar phi   = std::atan2(value[2], value[1]);

      value[1] = sigma;
      value[2] = tau;
      value[3] = phi;

      wrap_angles<target>(value);
    }
  }
  else if   constexpr (source == coordinate_system_type::cylindrical)
  {
    if      constexpr (target == coordinate_system_type::boyer_lindquist)
    {
      convert<source, coordinate_system_type::cartesian>(value);
      convert<coordinate_system_type::cartesian, target>(value, free_parameter);
    }
    else if constexpr (target == coordinate_system_type::prolate_spheroidal)
    {
      convert<source, coordinate_system_type::cartesian>(value);
      convert<coordinate_system_type::cartesian, target>(value, free_parameter);
    }
  }
  else if   constexpr (source == coordinate_system_type::prolate_spheroidal)
  {
    if      constexpr (target == coordinate_system_type::cartesian)
    {
      wrap_angles<source>(value);

      const scalar common = free_parameter * std::sqrt((static_cast<scalar>(std::pow(value[1], 2)) - static_cast<scalar>(1)) * (static_cast<scalar>(1) - static_cast<scalar>(std::pow(value[2], 2))));
      const scalar x      = common * std::cos(value[3]);
      const scalar y      = common * std::sin(value[3]);
      const scalar z      = free_parameter * value[1] * value[2];

      value[1] = x;
      value[2] = y;
      value[3] = z;
    }
    else if constexpr (target == coordinate_system_type::cylindrical)
    {
      convert<source, coordinate_system_type::cartesian>(value, free_parameter);
      convert<coordinate_system_type::cartesian, target>(value);
    }
    else if constexpr (target == coordinate_system_type::spherical)
    {
      convert<source, coordinate_system_type::cartesian>(value, free_parameter);
      convert<coordinate_system_type::cartesian, target>(value);
    }
  }
  else if   constexpr (source == coordinate_system_type::spherical)
  {
    if      constexpr (target == coordinate_system_type::boyer_lindquist)
    {
      convert<source, coordinate_system_type::cartesian>(value);
      convert<coordinate_system_type::cartesian, target>(value, free_parameter);
    }
    else if constexpr (target == coordinate_system_type::prolate_spheroidal)
    {
      convert<source, coordinate_system_type::cartesian>(value);
      convert<coordinate_system_type::cartesian, target>(value, free_parameter);
    }
  }
}

template<coordinate_system_type source, coordinate_system_type target, typename type, typename scalar = typename type::vector_type::value_type>
__device__ __host__ constexpr void convert_ray(type& value)
{
  if        constexpr (source == coordinate_system_type::cartesian)
  {
    if      constexpr (target == coordinate_system_type::cylindrical)
    {
      convert<source, target>(value.position);

      const scalar& r = value.position[1];
      const scalar& p = value.position[2];

      matrix44<scalar> transform;
      transform << 
        1,                 0,               0, 0,
        0,   std::cos(p)    , std::sin(p)    , 0,
        0, - std::sin(p) / r, std::cos(p) / r, 0,
        0,                 0,               0, 1;

      value.direction = transform * value.direction;
    }
    else if constexpr (target == coordinate_system_type::spherical)
    {
      convert<source, target>(value.position);
      
      const scalar& r = value.position[1];
      const scalar& t = value.position[2];
      const scalar& p = value.position[3];

      matrix44<scalar> transform;
      transform << 
        1,                               0,                             0,                0,
        0,   std::sin(t) * std::cos(p)    , std::sin(t) * std::sin(p)    ,  std::cos(t)    ,
        0,   std::cos(t) * std::cos(p) / r, std::cos(t) * std::sin(p) / r, -std::sin(t) / r,
        0, - std::sin(p) / std::sin(t) / r, std::cos(p) / std::sin(t) / r,                0;

      value.direction = transform * value.direction;
    }
  }
  else if   constexpr (source == coordinate_system_type::cylindrical)
  {
    if      constexpr (target == coordinate_system_type::cartesian)
    {
      const scalar& r = value.position[1];
      const scalar& p = value.position[2];

      matrix44<scalar> transform;
      transform << 
        1,           0,                0, 0,
        0, std::cos(p), -r * std::sin(p), 0,
        0, std::sin(p),  r * std::cos(p), 0,
        0,           0,                0, 1;

      convert<source, target>(value.position);

      value.direction = transform * value.direction;
    }
    else if constexpr (target == coordinate_system_type::spherical)
    {
      convert<source, target>(value.position);

      const scalar& r = value.position[1];
      const scalar& t = value.position[2];

      matrix44<scalar> transform;
      transform << 
        1,               0, 0,                0,
        0, std::sin(t)    , 0,  std::cos(t)    ,
        0, std::cos(t) / r, 0, -std::sin(t) / r,
        0,               0, 1,                0;

      value.direction = transform * value.direction;
    }
  }
  else if   constexpr (source == coordinate_system_type::spherical)
  {
    if      constexpr (target == coordinate_system_type::cartesian)
    {
      const scalar& r = value.position[1];
      const scalar& t = value.position[2];
      const scalar& p = value.position[3];

      matrix44<scalar> transform;
      transform << 
        1,                         0,                              0,                              0,
        0, std::sin(t) * std::cos(p),  r * std::cos(t) * std::cos(p), -r * std::sin(t) * std::sin(p),
        0, std::sin(t) * std::sin(p),  r * std::cos(t) * std::sin(p),  r * std::sin(t) * std::cos(p),
        0, std::cos(t)              , -r * std::sin(t)              ,                              0;

      convert<source, target>(value.position);

      value.direction = transform * value.direction;
    }
    else if constexpr (target == coordinate_system_type::cylindrical)
    {
      const scalar& r = value.position[1];
      const scalar& t = value.position[2];

      matrix44<scalar> transform;
      transform << 
        1,           0,                0, 0,
        0, std::sin(t),  r * std::cos(t), 0,
        0,           0,                0, 1,
        0, std::cos(t), -r * std::sin(t), 0;

      convert<source, target>(value.position);

      value.direction = transform * value.direction;
    }
  }
}
template<coordinate_system_type source, coordinate_system_type target, typename type, typename scalar = typename type::vector_type::value_type>
__device__ __host__ constexpr void convert_ray(type& value, const scalar free_parameter)
{
  if        constexpr (source == coordinate_system_type::boyer_lindquist)
  {
    if      constexpr (target == coordinate_system_type::cartesian)
    {
      const scalar& r = value.position[1];
      const scalar& t = value.position[2];
      const scalar& p = value.position[3];
      const scalar  k = std::sqrt(static_cast<scalar>(std::pow(free_parameter, 2)) + static_cast<scalar>(std::pow(r, 2)));

      matrix44<scalar> transform;
      transform << 
        1,                                 0,                             0,                              0,
        0, r * std::sin(t) * std::cos(p) / k, std::cos(t) * std::cos(p) * k, -std::sin(t) * std::sin(p) * k,
        0, r * std::sin(t) * std::sin(p) / k, std::cos(t) * std::sin(p) * k,  std::sin(t) * std::cos(p) * k,
        0,                       std::cos(t),              -r * std::sin(t),                              0;

      convert<source, target>(value.position, free_parameter);

      value.direction = transform * value.direction;
    }
    else if constexpr (target == coordinate_system_type::cylindrical)
    {
      convert<source, coordinate_system_type::cartesian>(value, free_parameter);
      convert<coordinate_system_type::cartesian, target>(value);
    }
    else if constexpr (target == coordinate_system_type::spherical)
    {
      convert<source, coordinate_system_type::cartesian>(value, free_parameter);
      convert<coordinate_system_type::cartesian, target>(value);
    }
  }
  else if   constexpr (source == coordinate_system_type::cartesian)
  {
    if      constexpr (target == coordinate_system_type::boyer_lindquist)
    {
      convert<source, target>(value.position, free_parameter);

      const scalar& r = value.position[1];
      const scalar& t = value.position[2];
      const scalar& p = value.position[3];
      const scalar  k = std::sqrt(static_cast<scalar>(std::pow(free_parameter, 2)) + static_cast<scalar>(std::pow(r, 2)));

      matrix44<scalar> transform;
      transform << 
        1,                                 0,                             0,                              0,
        0, r * std::sin(t) * std::cos(p) / k, std::cos(t) * std::cos(p) * k, -std::sin(t) * std::sin(p) * k,
        0, r * std::sin(t) * std::sin(p) / k, std::cos(t) * std::sin(p) * k,  std::sin(t) * std::cos(p) * k,
        0,                       std::cos(t),              -r * std::sin(t),                              0;
      transform = transform.inverse().eval(); // TODO: Ideally inversion should be analytic without runtime overhead. See https://github.com/einsteinpy/einsteinpy/blob/main/src/einsteinpy/coordinates/utils.py

      value.direction = transform * value.direction;
    }
    else if constexpr (target == coordinate_system_type::prolate_spheroidal)
    {
      convert<source, target>(value.position, free_parameter);

      const scalar& s = value.position[1];
      const scalar& t = value.position[2];
      const scalar& p = value.position[3];
      const scalar& a = free_parameter;
      const scalar  k = std::sqrt(-(static_cast<scalar>(std::pow(s, 2)) - static_cast<scalar>(1)) * (static_cast<scalar>(std::pow(t, 2)) - static_cast<scalar>(1)));

      matrix44<scalar> transform;
      transform << 
        1,                                                                    0,                                                                    0,                    0,
        0, -a * s * (static_cast<scalar>(std::pow(t, 2)) - 1) * std::cos(p) / k, -a * t * (static_cast<scalar>(std::pow(s, 2)) - 1) * std::cos(p) / k, -a * std::sin(p) * k,
        0, -a * s * (static_cast<scalar>(std::pow(t, 2)) - 1) * std::sin(p) / k, -a * t * (static_cast<scalar>(std::pow(s, 2)) - 1) * std::sin(p) / k,  a * std::cos(p) * k,
        0,  a * t                                                              ,  a * s                                                              ,                    0;
      transform = transform.inverse().eval(); // TODO: Ideally inversion should be analytic without runtime overhead.

      value.direction = transform * value.direction;
    }
  }
  else if   constexpr (source == coordinate_system_type::cylindrical)
  {
    if      constexpr (target == coordinate_system_type::boyer_lindquist)
    {
      convert<source, coordinate_system_type::cartesian>(value);
      convert<coordinate_system_type::cartesian, target>(value, free_parameter);
    }
    else if constexpr (target == coordinate_system_type::prolate_spheroidal)
    {
      convert<source, coordinate_system_type::cartesian>(value);
      convert<coordinate_system_type::cartesian, target>(value, free_parameter);
    }
  }
  else if   constexpr (source == coordinate_system_type::prolate_spheroidal)
  {
    if      constexpr (target == coordinate_system_type::cartesian)
    {
      const scalar& s = value.position[1];
      const scalar& t = value.position[2];
      const scalar& p = value.position[3];
      const scalar& a = free_parameter;
      const scalar  k = std::sqrt(-(static_cast<scalar>(std::pow(s, 2)) - 1) * (static_cast<scalar>(std::pow(t, 2)) - 1));

      matrix44<scalar> transform;
      transform << 
        1,                                                                    0,                                                                    0,                    0,
        0, -a * s * (static_cast<scalar>(std::pow(t, 2)) - 1) * std::cos(p) / k, -a * t * (static_cast<scalar>(std::pow(s, 2)) - 1) * std::cos(p) / k, -a * std::sin(p) * k,
        0, -a * s * (static_cast<scalar>(std::pow(t, 2)) - 1) * std::sin(p) / k, -a * t * (static_cast<scalar>(std::pow(s, 2)) - 1) * std::sin(p) / k,  a * std::cos(p) * k,
        0,  a * t                                                              ,  a * s                                                              ,                    0;

      convert<source, target>(value.position, free_parameter);

      value.direction = transform * value.direction;
    }
    else if constexpr (target == coordinate_system_type::cylindrical)
    {
      convert<source, coordinate_system_type::cartesian>(value, free_parameter);
      convert<coordinate_system_type::cartesian, target>(value);
    }
    else if constexpr (target == coordinate_system_type::spherical)
    {
      convert<source, coordinate_system_type::cartesian>(value, free_parameter);
      convert<coordinate_system_type::cartesian, target>(value);
    }
  }
  else if   constexpr (source == coordinate_system_type::spherical)
  {
    if      constexpr (target == coordinate_system_type::boyer_lindquist)
    {
      convert<source, coordinate_system_type::cartesian>(value);
      convert<coordinate_system_type::cartesian, target>(value, free_parameter);
    }
    else if constexpr (target == coordinate_system_type::prolate_spheroidal)
    {
      convert<source, coordinate_system_type::cartesian>(value);
      convert<coordinate_system_type::cartesian, target>(value, free_parameter);
    }
  }
}
}