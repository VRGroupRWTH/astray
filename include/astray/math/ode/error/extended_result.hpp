#pragma once

namespace ast
{
template <typename type_>
struct extended_result
{
  using type = type_;

  type value;
  type error;
};
}