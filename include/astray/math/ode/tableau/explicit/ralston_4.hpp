#pragma once

#include <astray/math/ode/tableau/tableau_traits.hpp>

namespace ast
{
template <typename type = double>
struct ralston_4_tableau {};

template <typename type>
constexpr __constant__ auto tableau_a<ralston_4_tableau<type>> = std::array
{
  type(0.4),
  type(0.29697761), type( 0.15875964),
  type(0.21810040), type(-3.05096516), type(3.83286476)
};
template <typename type>
constexpr __constant__ auto tableau_b<ralston_4_tableau<type>> = std::array
{
  type(0.17476028), type(-0.55148066), type(1.20553560), type(0.17118478)
};
template <typename type>
constexpr __constant__ auto tableau_c<ralston_4_tableau<type>> = std::array
{
  type(0.0), type(0.4), type(0.45573725), type(1.0)
};
}