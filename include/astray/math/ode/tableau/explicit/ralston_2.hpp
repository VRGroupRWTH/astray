#pragma once

#include <astray/math/ode/tableau/tableau_traits.hpp>

namespace ast
{
template <typename type = double>
struct ralston_2_tableau {};

template <typename type>
constexpr __constant__ auto tableau_a<ralston_2_tableau<type>> = std::array {type(2.0 / 3.0)};
template <typename type>
constexpr __constant__ auto tableau_b<ralston_2_tableau<type>> = std::array {type(0.25), type(0.75)};
template <typename type>
constexpr __constant__ auto tableau_c<ralston_2_tableau<type>> = std::array {type(0.0) , type(2.0 / 3.0)};
}