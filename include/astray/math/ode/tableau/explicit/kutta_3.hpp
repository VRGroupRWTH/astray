#pragma once

#include <astray/math/ode/tableau/tableau_traits.hpp>

namespace ast
{
template <typename type = double>
struct kutta_3_tableau {};

template <typename type>
constexpr __constant__ auto tableau_a<kutta_3_tableau<type>> = std::array
{
  type( 0.5),
  type(-1.0), type(2.0)
};
template <typename type>
constexpr __constant__ auto tableau_b<kutta_3_tableau<type>> = std::array {type(1.0 / 6.0), type(2.0 / 3.0), type(1.0 / 6.0)};
template <typename type>
constexpr __constant__ auto tableau_c<kutta_3_tableau<type>> = std::array {type(0.0), type(0.5), type(1.0)};
}