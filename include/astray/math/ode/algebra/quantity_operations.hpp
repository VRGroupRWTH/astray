#pragma once

#include <cstddef>

#include <astray/math/ode/algebra/quantity_traits.hpp>

namespace ast
{
template <typename value_type>
struct quantity_operations 
{
  template <typename function_type, typename... value_types>
  static constexpr void        for_each(const function_type& function, const value_types&... values)
  {
    if constexpr (is_scalar_v<value_type>)
      function(values...);
    else
    {
      const auto element_count = std::get<0>(std::tie(values...)).size();
      for (auto i = 0; i < element_count; ++i)
        function(values.data()[i]...);
    }
  }
  static constexpr std::size_t size    (const value_type& value)
  {
    if constexpr (is_scalar_v<value_type>)
      return 1;
    else
      return value.size();
  }
};
}