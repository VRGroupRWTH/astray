#pragma once

#include <complex>
#include <cstdint>
#include <type_traits>

#include <astray/math/linear_algebra.hpp>

namespace ast
{
template <typename type>
inline constexpr bool is_complex_v                                                                    = false;
template <typename type>
inline constexpr bool is_complex_v<std::complex<type>>                                                = true ;

template <typename type>
inline constexpr bool is_scalar_v                                                                     = std::is_arithmetic_v<type> || is_complex_v<type>;

template <typename type>
inline constexpr bool is_vector_v                                                                     = false;
template <typename type, std::int32_t size, std::int32_t options, std::int32_t max_rows>
inline constexpr bool is_vector_v<Eigen::Matrix<type, size, 1, options, max_rows, 1>>                 = true ;
template <typename type, std::int32_t size, std::int32_t options, std::int32_t max_columns>
inline constexpr bool is_vector_v<Eigen::Matrix<type, 1, size, options, 1, max_columns>>              = true ;

template <typename type>
inline constexpr bool is_matrix_v                                                                     = false;
template <typename type, std::int32_t rows, std::int32_t columns, std::int32_t options, std::int32_t max_rows, std::int32_t max_columns>
inline constexpr bool is_matrix_v<Eigen::Matrix<type, rows, columns, options, max_rows, max_columns>> = rows != 1 && columns != 1;

template <typename type>
inline constexpr bool is_tensor_v                                                                     = false;
template <typename type, std::int32_t sizes, std::int32_t options, typename index_type>
inline constexpr bool is_tensor_v<Eigen::Tensor         <type, sizes, options, index_type>>           = true ;
template <typename type, typename     sizes, std::int32_t options, typename index_type>
inline constexpr bool is_tensor_v<Eigen::TensorFixedSize<type, sizes, options, index_type>>           = true ;

template <typename type>
struct is_complex : std::bool_constant<is_complex_v<type>> { };
template <typename type>
struct is_scalar  : std::bool_constant<is_scalar_v <type>> { };
template <typename type>
struct is_vector  : std::bool_constant<is_vector_v <type>> { };
template <typename type>
struct is_matrix  : std::bool_constant<is_matrix_v <type>> { };
template <typename type>
struct is_tensor  : std::bool_constant<is_tensor_v <type>> { };
}