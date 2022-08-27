#pragma once

#include <cstddef>

namespace ast
{
template <typename type, std::size_t unused>
using parameter_pack_expander = type;
}