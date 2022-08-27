#pragma once

#include <variant>

#include <astray/math/angle.hpp>

namespace ast
{
template <typename type>
struct perspective_projection
{
  type fov_y        = to_radians <type>(60  );
  type aspect_ratio = static_cast<type>(1   );
  type focal_length = static_cast<type>(1   );
  type near_clip    = static_cast<type>(0.01);
  type far_clip     = static_cast<type>(100 );
};

template <typename type>
struct orthographic_projection
{
  type height       = static_cast<type>(1   );
  type aspect_ratio = static_cast<type>(1   );
  type near_clip    = static_cast<type>(0.01);
  type far_clip     = static_cast<type>(100 );
};

template <typename type>
using projection = std::variant<
  perspective_projection <type>, 
  orthographic_projection<type>>;
}