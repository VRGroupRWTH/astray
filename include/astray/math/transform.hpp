#pragma once

#include <astray/math/linear_algebra.hpp>

namespace ast
{
template <typename scalar_type_>
struct transform
{
  using scalar_type     = scalar_type_;
  using vector_type     = vector3   <scalar_type>;
  using matrix_type     = matrix33  <scalar_type>;
  using quaternion_type = quaternion<scalar_type>;
  using angle_axis_type = angle_axis<scalar_type>;

  vector_type right              () const
  {
    return rotation * vector_type::UnitX();
  }
  vector_type up                 () const
  {
    return rotation * vector_type::UnitY();
  }
  vector_type forward            () const
  {
    return rotation * vector_type::UnitZ();
  }

  vector_type rotation_to_euler  () const
  {
    return rotation.toRotationMatrix().eulerAngles(0, 1, 2);
  }
  void        rotation_from_euler(const vector_type& euler)
  {
    rotation = angle_axis_type(euler[0], vector_type::UnitX())
             * angle_axis_type(euler[1], vector_type::UnitY())
             * angle_axis_type(euler[2], vector_type::UnitZ());
  }

  void        look_at            (const vector_type& target, const vector_type& up = vector_type::UnitY())
  {
    const auto forward = (target - translation).normalized();
    const auto right   = up     .cross(forward).normalized();
    const auto true_up = forward.cross(right  ).normalized();

    matrix_type matrix;
    matrix << right, true_up, forward;

    rotation = quaternion_type(matrix);
  }
  void        clear              ()
  {
    translation.setZero    ();
    rotation   .setIdentity();
    scale      .setOnes    ();
  }

  vector_type     translation = vector_type    ::Zero    ();
  quaternion_type rotation    = quaternion_type::Identity();
  vector_type     scale       = vector_type    ::Ones    ();
};
}