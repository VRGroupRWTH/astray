#pragma once

#define EIGEN_INITIALIZE_MATRICES_BY_ZERO

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <unsupported/Eigen/CXX11/Tensor>

namespace ast
{
// Vectors.
template <typename type, std::int32_t size>
using vector     = Eigen::Matrix<type, size, 1>;
template <typename type>
using vector2    = vector<type, 2>;
template <typename type>
using vector3    = vector<type, 3>;
template <typename type>
using vector4    = vector<type, 4>;

// Matrices.
template <typename type, std::int32_t rows, std::int32_t columns>
using matrix     = Eigen::Matrix<type, rows, columns>;
template <typename type>
using matrix22   = matrix<type, 2, 2>;
template <typename type>
using matrix33   = matrix<type, 3, 3>;
template <typename type>
using matrix44   = matrix<type, 4, 4>;

// Tensors.
template <typename type, typename sizes>
using tensor     = Eigen::TensorFixedSize<type, sizes>;
template <typename type>
using tensor222  = tensor<type, Eigen::Sizes<2, 2, 2>>;
template <typename type>
using tensor333  = tensor<type, Eigen::Sizes<3, 3, 3>>;
template <typename type>
using tensor444  = tensor<type, Eigen::Sizes<4, 4, 4>>;

// AABBs.
template <typename type, std::int32_t size>
using aabb       = Eigen::AlignedBox<type, size>;
template <typename type>
using aabb2      = aabb<type, 2>;
template <typename type>
using aabb3      = aabb<type, 3>;
template <typename type>
using aabb4      = aabb<type, 4>;

// Quaternions.
template <typename type>
using quaternion = Eigen::Quaternion<type>;

// Angle axes.
template <typename type>
using angle_axis = Eigen::AngleAxis <type>;

// Mapped types.
template <typename type>
using mapped     = Eigen::Map       <type>;
}