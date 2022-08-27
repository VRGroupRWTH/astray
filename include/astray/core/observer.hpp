#pragma once

#include <astray/math/indexing.hpp>
#include <astray/math/projection.hpp>
#include <astray/math/ray.hpp>
#include <astray/math/transform.hpp>
#include <astray/parallel/thrust.hpp>

namespace ast
{
template <typename scalar_type>
class observer
{
public:
  using transform_type  = transform <scalar_type>;
  using projection_type = projection<scalar_type>;
  using ray_type        = ray<vector4<scalar_type>>;

  using vector_type     = typename transform_type::vector_type;
  using image_size_type = vector2<std::int32_t>;

  struct device_data_perspective
  {
    vector_type     position       ;
    vector_type     direction_00   ;
    vector_type     u              ;
    vector_type     v              ;
    image_size_type global_size    ;
    image_size_type local_size     ;
    image_size_type local_offset   ;
    scalar_type     coordinate_time;
  };
  struct device_data_orthographic
  {
    vector_type     direction      ;
    vector_type     position_00    ;
    vector_type     u              ;
    vector_type     v              ;
    image_size_type global_size    ;
    image_size_type local_size     ;
    image_size_type local_offset   ;
    scalar_type     coordinate_time;
  };
  
  thrust::device_vector<ray_type>& generate_rays(
    const image_size_type& global_size ,
    const image_size_type& local_size  ,
    const image_size_type& local_offset)
  {
    const auto linear_size = local_size.prod();
    if (rays_.size() != linear_size)
      rays_.resize(linear_size);

    if (std::holds_alternative<perspective_projection<scalar_type>>(projection_))
      return generate_rays_perspective (global_size, local_size, local_offset);
    else
      return generate_rays_orthographic(global_size, local_size, local_offset);
  }
  
  thrust::device_vector<ray_type>& generate_rays_perspective(
    const image_size_type& global_size ,
    const image_size_type& local_size  ,
    const image_size_type& local_offset)
  {
    const auto&       cast_projection = std::get<perspective_projection<scalar_type>>(projection_);
    const scalar_type v_size          = static_cast<scalar_type>(2) * tan(static_cast<scalar_type>(0.5) * cast_projection.fov_y);
    const scalar_type u_size          = v_size * cast_projection.aspect_ratio;

    const vector_type u               = transform_.right  () * u_size;
    const vector_type v               = transform_.up     () * v_size;
    const vector_type w               = transform_.forward() * cast_projection.focal_length;
    const vector_type direction_00    = w - static_cast<scalar_type>(0.5) * u + static_cast<scalar_type>(0.5) * v;

    device_data_perspective device_data
    {
      transform_.translation,
      direction_00         ,
      u                    ,
      v                    ,
      global_size          ,
      local_size           ,
      local_offset         ,
      coordinate_time_
    };
    thrust::copy_n(&device_data, 1, perspective_data_.begin());

    thrust::for_each(
      thrust::make_zip_iterator(thrust::make_tuple(thrust::counting_iterator<std::size_t>(0)           , rays_.begin())),
      thrust::make_zip_iterator(thrust::make_tuple(thrust::counting_iterator<std::size_t>(rays_.size()), rays_.end  ())),
      [data = perspective_data_.data().get()] __device__ (const auto& iteratee)
      {
        const auto local_index        = thrust::get<0>(iteratee);
        const auto local_multi_index  = unravel_index<image_size_type, true>(local_index, data->local_size);
        const auto global_multi_index = local_multi_index + data->local_offset;

        const auto direction          = (data->direction_00
          + data->u * static_cast<scalar_type>(global_multi_index[0]) / static_cast<scalar_type>(data->global_size[0])
          - data->v * static_cast<scalar_type>(global_multi_index[1]) / static_cast<scalar_type>(data->global_size[1])).normalized();

        auto& ray             = thrust::get<1>(iteratee);
        ray.position [0]      = data->coordinate_time;
        ray.position .tail(3) = data->position;
        ray.direction[0]      = static_cast<scalar_type>(-1);
        ray.direction.tail(3) = direction;
      });

    return rays_;
  }

  thrust::device_vector<ray_type>& generate_rays_orthographic(
    const image_size_type& global_size ,
    const image_size_type& local_size  ,
    const image_size_type& local_offset)
  {
    const auto&       cast_projection = std::get<orthographic_projection<scalar_type>>(projection_);
    const vector_type u               = transform_.right() * cast_projection.height * cast_projection.aspect_ratio;
    const vector_type v               = transform_.up   () * cast_projection.height;
    const vector_type position_00     = transform_.translation - static_cast<scalar_type>(0.5) * u + static_cast<scalar_type>(0.5) * v;

    device_data_orthographic device_data
    {
      transform_.forward(),
      position_00        ,
      u                  ,
      v                  ,
      global_size        ,
      local_size         ,
      local_offset       ,
      coordinate_time_
    };
    thrust::copy_n(&device_data, 1, orthographic_data_.begin());
    
    thrust::for_each(
      thrust::make_zip_iterator(thrust::make_tuple(thrust::counting_iterator<std::size_t>(0)           , rays_.begin())),
      thrust::make_zip_iterator(thrust::make_tuple(thrust::counting_iterator<std::size_t>(rays_.size()), rays_.end  ())),
      [data = orthographic_data_.data().get()] __device__ (const auto& iteratee)
      {
        const auto local_index        = thrust::get<0>(iteratee);
        const auto local_multi_index  = unravel_index<image_size_type, true>(local_index, data->local_size);
        const auto global_multi_index = local_multi_index + data->local_offset;
        
        const auto position           = data->position_00
          + data->u * static_cast<scalar_type>(global_multi_index[0]) / static_cast<scalar_type>(data->global_size[0])
          - data->v * static_cast<scalar_type>(global_multi_index[1]) / static_cast<scalar_type>(data->global_size[1]);

        auto& ray             = thrust::get<1>(iteratee);
        ray.position [0]      = data->coordinate_time;
        ray.position .tail(3) = position;
        ray.direction[0]      = static_cast<scalar_type>(-1);
        ray.direction.tail(3) = data->direction;
      });

    return rays_;
  }
  
  scalar_type                      get_coordinate_time() const
  {
    return coordinate_time_;
  }
  void                             set_coordinate_time(const scalar_type      value)
  {
    coordinate_time_ = value;
  }
  
        transform_type&            get_transform      ()
  {
    return transform_;
  }
  const transform_type&            get_transform      () const
  {
    return transform_;
  }
  void                             set_transform      (const transform_type&  value)
  {
    transform_ = value;
  }

        projection_type&           get_projection     ()
  {
    return projection_;
  }
  const projection_type&           get_projection     () const
  {
    return projection_;
  }
  void                             set_projection     (const projection_type& value)
  {
    projection_ = value;
  }

protected:
  scalar_type                                     coordinate_time_   = static_cast<scalar_type>(0);
  transform_type                                  transform_         { };
  projection_type                                 projection_        { };

  thrust::device_vector<ray_type>                 rays_              { };
  thrust::device_vector<device_data_perspective>  perspective_data_  {1};
  thrust::device_vector<device_data_orthographic> orthographic_data_ {1};
};
}