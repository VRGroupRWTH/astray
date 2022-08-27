#pragma once

#include <astray/math/ode/tableau/tableau_traits.hpp>

namespace ast
{
template <typename type = double>
struct heun_2_tableau {};

template <typename type>
constexpr __constant__ auto tableau_a<heun_2_tableau<type>> = std::array {type(1.0)};
template <typename type>
constexpr __constant__ auto tableau_b<heun_2_tableau<type>> = std::array {type(0.5), type(0.5)};
template <typename type>
constexpr __constant__ auto tableau_c<heun_2_tableau<type>> = std::array {type(0.0), type(1.0)};
}