#pragma once

#include <astray/math/ode/tableau/tableau_traits.hpp>

namespace ast
{
template <typename type = double>
struct ralston_3_tableau {};

template <typename type>
constexpr __constant__ auto tableau_a<ralston_3_tableau<type>> = std::array
{
  type(0.5),
  type(0.0), type(0.75)
};
template <typename type>
constexpr __constant__ auto tableau_b<ralston_3_tableau<type>> = std::array {type(2.0 / 9.0), type(1.0 / 3.0), type(4.0 / 9.0)};
template <typename type>
constexpr __constant__ auto tableau_c<ralston_3_tableau<type>> = std::array {type(0.0), type(0.5), type(0.75)};
}