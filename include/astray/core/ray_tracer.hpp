#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include <astray/core/geodesic.hpp>
#include <astray/core/observer.hpp>
#include <astray/math/constants.hpp>
#include <astray/math/coordinate_system.hpp>
#include <astray/media/image.hpp>
#include <astray/parallel/mpi/mpi.hpp>
#include <astray/parallel/partitioner.hpp>
#include <astray/parallel/thrust.hpp>

namespace ast
{
template <typename metric_type, typename motion_type>
class ray_tracer
{
public:
  using scalar_type          = typename motion_type::scalar_type;
  using vector_type          = vector4<scalar_type>;

  using observer_type        = observer<scalar_type>;

  using pixel_type           = vector3<std::uint8_t>;
  using image_type           = image<pixel_type>;
  using image_size_type      = image_type::size_type;

  using bounds_type          = typename motion_type::bounds_type;
  using error_evaluator_type = typename motion_type::error_evaluator_type;

  using partitioner_type     = partitioner<2, std::int32_t, image_size_type, true>;

  struct device_data
  {
    vector_type          observer_position ;
    pixel_type*          background        ;
    image_size_type      background_size   ;

    metric_type          metric            ;
    std::size_t          iterations        ;
    scalar_type          lambda_step_size  ;
    scalar_type          lambda            ;
    bounds_type          bounds            ;
    //error_evaluator_type error_evaluator   ;
    bool                 debug             ;
    
    pixel_type*          result            ;
    image_size_type      result_size       ;
    image_size_type      result_offset     ;
  };

  explicit ray_tracer  (
    const image_size_type&      image_size       = {1920, 1080},
    const metric_type&          metric           = metric_type(),
    const std::size_t           iterations       = static_cast<std::size_t>(1e3),
    const scalar_type           lambda_step_size = static_cast<scalar_type>(1e-3),
    const scalar_type           lambda           = static_cast<scalar_type>(0),
    const bounds_type&          bounds           = bounds_type(),
    const error_evaluator_type& error_evaluator  = error_evaluator_type(),
    const bool                  debug            = false)
  : metric_            (metric)
  , iterations_        (iterations)
  , lambda_step_size_  (lambda_step_size)
  , lambda_            (lambda)
  , bounds_            (bounds)
  , error_evaluator_   (error_evaluator)
  , debug_             (debug)
  , device_background_ (background_.data)
  , partitioner_       (communicator_.rank(), communicator_.size(), image_size)
  {
#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
    constexpr auto target_heap_size = static_cast<std::size_t>(1e+9);
    std::size_t    heap_size        = 0;
    cudaDeviceGetLimit(&heap_size, cudaLimitMallocHeapSize);
    if (heap_size < target_heap_size)
      cudaDeviceSetLimit(cudaLimitMallocHeapSize, target_heap_size);
#endif

#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_OMP
    omp_set_num_threads(omp_get_num_procs());
#endif

    set_image_size(image_size);
  }
  ray_tracer           (const ray_tracer&  that) = delete ;
  ray_tracer           (      ray_tracer&& temp) = default;
 ~ray_tracer           ()                        = default;
  ray_tracer& operator=(const ray_tracer&  that) = delete ;
  ray_tracer& operator=(      ray_tracer&& temp) = default;
  
  const image_type&           render_frame        ()
  {
    using constants = constants<scalar_type>;

    device_data data 
    {
      vector_type(
        observer_.get_coordinate_time(), 
        observer_.get_transform().translation[0], 
        observer_.get_transform().translation[1], 
        observer_.get_transform().translation[2]),
      device_background_.data().get(),
      background_.size               ,
      metric_                        ,
      iterations_                    ,
      lambda_step_size_              ,
      lambda_                        ,
      bounds_                        ,
      //error_evaluator              ,
      debug_                         ,
      device_result_.data().get()    ,
      result_.size                   ,
      partitioner_.rank_offset()
    };
    thrust::copy_n(&data, 1, device_data_.begin());
    
    auto& rays = observer_.generate_rays(partitioner_.domain_size(), partitioner_.block_size(), partitioner_.rank_offset());

    thrust::for_each(
      thrust::make_zip_iterator(thrust::make_tuple(thrust::counting_iterator<std::size_t>(0)          , rays.begin())),
      thrust::make_zip_iterator(thrust::make_tuple(thrust::counting_iterator<std::size_t>(rays.size()), rays.end  ())),
      [data = device_data_.data().get()] __device__ (const auto& iteratee)
      {
        auto  index = thrust::get<0>(iteratee);
        auto& ray   = thrust::get<1>(iteratee);
        
        metric_type metric(data->metric); // A copy is necessary for correct creation of virtual function table.
        
        if constexpr (metric_type::coordinate_system() == coordinate_system_type::boyer_lindquist || metric_type::coordinate_system() == coordinate_system_type::prolate_spheroidal)
          convert_ray<coordinate_system_type::cartesian, metric_type::coordinate_system()>(ray, metric.coordinate_system_parameter());
        else
          convert_ray<coordinate_system_type::cartesian, metric_type::coordinate_system()>(ray);
        
        const auto termination = motion_type::integrate(ray, metric, data->iterations, data->lambda_step_size, data->lambda, data->bounds); //, data->error_evaluator);
        
        if (termination == termination_reason::none || termination == termination_reason::out_of_bounds)
        {
          if constexpr (metric_type::coordinate_system() == coordinate_system_type::boyer_lindquist || metric_type::coordinate_system() == coordinate_system_type::prolate_spheroidal)
            convert<metric_type::coordinate_system(), coordinate_system_type::cartesian>(ray.position, metric.coordinate_system_parameter());
          else
            convert<metric_type::coordinate_system(), coordinate_system_type::cartesian>(ray.position);
        
          ray.position -= data->observer_position; // Environment map is relative to observer.
          
          convert<coordinate_system_type::cartesian, coordinate_system_type::spherical>(ray.position);
          
          image_size_type background_index(
            std::floor(ray.position[3] / constants::two_pi * static_cast<scalar_type>(data->background_size[0])),
            std::floor(ray.position[2] / constants::pi     * static_cast<scalar_type>(data->background_size[1])));
          if (background_index[1] == data->background_size[1]) --background_index[1];
          
          data->result[index] = data->background[ravel_multi_index<image_size_type, true>(background_index, data->background_size)];
        }
        
        if (data->debug)
        {
          if      (termination == termination_reason::constraint_violation)
            data->result[index] = pixel_type(255, 128, 128);
          else if (termination == termination_reason::numeric_error       )
            data->result[index] = pixel_type(128, 255, 128);
          else if (termination == termination_reason::spacetime_breakdown )
            data->result[index] = pixel_type(128, 128, 255);
        }
      });

    thrust::copy(device_result_.begin(), device_result_.end(), result_.data.begin());

#ifdef ASTRAY_USE_MPI
    std::vector<std::int32_t> counts         (communicator_.size(), 1);
    std::vector<std::int32_t> displacements  (communicator_.size());
    for (auto y = 0; y < partitioner_.grid_size()[1]; ++y)
      for (auto x = 0; x < partitioner_.grid_size()[0]; ++x)
        displacements[x + y * partitioner_.grid_size()[0]] = x + y * (partitioner_.block_size()[1] * partitioner_.grid_size()[0]);

    communicator_.gatherv(
      result_         .data.data(), static_cast<std::int32_t>(result_.data.size()), pixel_data_type_  ,
      gathered_result_.data.data(), counts.data(), displacements.data()           , resized_data_type_);

    if (communicator_.rank() != 0)
      return result_; // Workers return their partial results.
    return gathered_result_;
#else
    return result_;
#endif
  }

  const image_size_type&      get_image_size      () const
  {
    return partitioner_.domain_size();
  }
  void                        set_image_size      (const image_size_type&      value)
  {
    partitioner_.set_domain_size(value);
    
    result_             = image_type(partitioner_.block_size(), pixel_type());
    device_result_      = result_.data;

#ifdef ASTRAY_USE_MPI
    gathered_result_    = image_type(partitioner_.domain_size());

    pixel_data_type_    = mpi::data_type(mpi::data_type(MPI_UNSIGNED_CHAR), 3);
    subarray_data_type_ = mpi::data_type(pixel_data_type_   , value, partitioner_.block_size(), image_size_type::Zero().eval(), true);
    resized_data_type_  = mpi::data_type(subarray_data_type_, 0    , partitioner_.block_size()[0] * sizeof(pixel_type));
#endif
    
    std::visit([value] (auto& visited) 
    { 
      visited.aspect_ratio = static_cast<scalar_type>(value[0]) / static_cast<scalar_type>(value[1]); 
    }, observer_.get_projection());
  }

        observer_type&        get_observer        ()
  {
    return observer_;
  }
  const observer_type&        get_observer        () const
  {
    return observer_;
  }
  void                        set_observer        (const observer_type&        value)
  {
    observer_ = value;
  }

  const image_type&           get_background      () const
  {
    return background_;
  }
  void                        set_background      (const image_type&           value)
  {
    background_        = value;
    device_background_ = background_.data;
  }

        metric_type&          get_metric          ()
  {
    return metric_;
  }
  const metric_type&          get_metric          () const
  {
    return metric_;
  }
  void                        set_metric          (const metric_type&          value)
  {
    metric_ = value;
  }

  std::size_t                 get_iterations      () const
  {
    return iterations_;
  }
  void                        set_iterations      (const std::size_t           value)
  {
    iterations_ = value;
  }

  scalar_type                 get_lambda_step_size() const
  {
    return lambda_step_size_;
  }
  void                        set_lambda_step_size(const scalar_type           value)
  {
    lambda_step_size_ = value;
  }

  scalar_type                 get_lambda          () const
  {
    return lambda_;
  }
  void                        set_lambda          (const scalar_type           value)
  {
    lambda_ = value;
  }

        bounds_type&          get_bounds          ()
  {
    return bounds_;
  }
  const bounds_type&          get_bounds          () const
  {
    return bounds_;
  }
  void                        set_bounds          (const bounds_type&          value)
  {
    bounds_ = value;
  }

        error_evaluator_type& get_error_evaluator ()
  {
    return error_evaluator_;
  }
  const error_evaluator_type& get_error_evaluator () const
  {
    return error_evaluator_;
  }
  void                        set_error_evaluator (const error_evaluator_type& value)
  {
    error_evaluator_ = value;
  }
  
  bool                        is_debug            () const
  {
    return debug_;
  }
  void                        set_debug           (const bool                  value)
  {
    debug_ = value;
  }

  const mpi::communicator&    get_communicator    () const
  {
    return communicator_;
  }
  
protected:
  observer_type                      observer_          ;
  image_type                         background_        ;
                                     
  metric_type                        metric_            ;
  std::size_t                        iterations_        ;
  scalar_type                        lambda_step_size_  ;
  scalar_type                        lambda_            ;
  bounds_type                        bounds_            ;
  error_evaluator_type               error_evaluator_   ;
  bool                               debug_             ;

  thrust::device_vector<device_data> device_data_       {1};
  thrust::device_vector<pixel_type>  device_background_ ;
  thrust::device_vector<pixel_type>  device_result_     ;
  image_type                         result_            ;
  image_type                         gathered_result_   ;

  mpi::environment                   environment_       ;
  mpi::communicator                  communicator_      ;
  partitioner_type                   partitioner_       ;
  mpi::data_type                     pixel_data_type_   ;
  mpi::data_type                     subarray_data_type_;
  mpi::data_type                     resized_data_type_ ;
};
}