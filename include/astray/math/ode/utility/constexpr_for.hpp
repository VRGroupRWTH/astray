#pragma once

#include <type_traits>

namespace ast
{
template <auto begin, auto end, auto increment, typename function_type>
constexpr void constexpr_for(function_type&& function)
{
  if constexpr (begin < end)
  {
    function(std::integral_constant<decltype(begin), begin>());
    constexpr_for<begin + increment, end, increment>(function);
  }
}
}