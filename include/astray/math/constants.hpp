#pragma once

#define _USE_MATH_DEFINES

#include <cmath>
#include <limits>
#include <math.h>

namespace ast
{
template <typename type>
struct constants
{
static constexpr type speed_of_light             = static_cast<type>(1);
static constexpr type speed_of_light_squared     = static_cast<type>(speed_of_light * speed_of_light);

static constexpr type gravitational_constant     = static_cast<type>(1);
static constexpr type two_gravitational_constant = static_cast<type>(2) * gravitational_constant;

static constexpr type vacuum_permeability        = static_cast<type>(1);
static constexpr type vacuum_permittivity        = static_cast<type>(1) / (vacuum_permeability * speed_of_light_squared);

static constexpr type cosmological_constant      = static_cast<type>(1);

static constexpr type pi                         = static_cast<type>(M_PI);
static constexpr type two_pi                     = static_cast<type>(2) * pi;
static constexpr type three_pi                   = static_cast<type>(3) * pi;
static constexpr type four_pi                    = static_cast<type>(4) * pi;

static constexpr type epsilon                    = std::numeric_limits<type>::epsilon();

static constexpr type schwarzschild_radius       (const type mass  )
{
  return two_gravitational_constant * mass / speed_of_light_squared;
}
static constexpr type characteristic_length_scale(const type charge)
{
  return std::pow(charge, 2) * gravitational_constant / (vacuum_permittivity * std::pow(speed_of_light, 4));
}
};
}