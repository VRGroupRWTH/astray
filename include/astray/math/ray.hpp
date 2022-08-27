#pragma once

#include <astray/math/linear_algebra.hpp>

namespace ast
{
template <typename vector_type_>
struct ray
{
  using vector_type = vector_type_;

  vector_type position  = vector_type::Zero();
  vector_type direction = vector_type::Zero();
};
}