#pragma once

#include <astray/math/ode/tableau/tableau_traits.hpp>

namespace ast
{
template <typename type = double>
struct forward_euler_tableau {};

template <typename type>
constexpr __constant__ auto tableau_a<forward_euler_tableau<type>> = std::array {type(0.0)};
template <typename type>
constexpr __constant__ auto tableau_b<forward_euler_tableau<type>> = std::array {type(1.0)};
template <typename type>
constexpr __constant__ auto tableau_c<forward_euler_tableau<type>> = std::array {type(0.0)};
}