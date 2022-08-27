#pragma once

#include <cstdint>
#include <vector>

#ifdef ASTRAY_USE_MPI
#include <mpi.h>
#endif

namespace ast::mpi
{
class data_type
{
public:
  explicit data_type  (const std::int32_t native = 0, const bool managed = false) 
  : native_(native), managed_(managed)
  {
  
  }
  data_type           (const data_type&  that, const std::int32_t size) 
  : managed_(true)
  {
#ifdef ASTRAY_USE_MPI
    MPI_Type_contiguous(size, that.native(), &native_);
    MPI_Type_commit    (&native_);
#endif
  }
  data_type           (const data_type&  that, const std::int64_t lower_bound, const std::int64_t extent) 
  : managed_(true)
  {
#ifdef ASTRAY_USE_MPI
    MPI_Type_create_resized(that.native(), lower_bound, extent, &native_);
    MPI_Type_commit        (&native_);
#endif
  }
  template <typename type>
  data_type           (const data_type&  that, const type& sizes, const type& sub_sizes, const type& starts, const bool fortran_order = false)
  : managed_(true)
  {
#ifdef ASTRAY_USE_MPI
    MPI_Type_create_subarray(static_cast<std::int32_t>(sizes.size()), sizes.data(), sub_sizes.data(), starts.data(), fortran_order ? MPI_ORDER_FORTRAN : MPI_ORDER_C, that.native(), &native_);
    MPI_Type_commit         (&native_);
#endif
  }
  data_type           (const data_type&  that)
  : managed_(true)
  {
#ifdef ASTRAY_USE_MPI
    MPI_Type_dup(that.native_, &native_);
#endif
  }
  data_type           (      data_type&& temp) noexcept
  : native_(temp.native_), managed_(temp.managed_)
  {
#ifdef ASTRAY_USE_MPI
    temp.managed_ = false;
    temp.native_  = MPI_DATATYPE_NULL;
#endif
  }
  virtual ~data_type  ()
  {
#ifdef ASTRAY_USE_MPI
    if (managed_)
      MPI_Type_free(&native_);
#endif
  }
  data_type& operator=(const data_type&  that)
  {
    if (this != &that)
    {
#ifdef ASTRAY_USE_MPI
      if (managed_ && native_ != MPI_DATATYPE_NULL)
        MPI_Type_free(&native_);

      managed_ = true;
      MPI_Type_dup(that.native_, &native_);
#endif
    }
    return *this;
  }
  data_type& operator=(      data_type&& temp) noexcept
  {
    if (this != &temp)
    {
#ifdef ASTRAY_USE_MPI
      if (managed_ && native_ != MPI_DATATYPE_NULL)
        MPI_Type_free(&native_);
      
      managed_      = temp.managed_;
      native_       = temp.native_ ;

      temp.managed_ = false;
      temp.native_  = MPI_DATATYPE_NULL;
#endif
    }
    return *this;
  }
  
  constexpr std::int32_t native() const
  {
    return native_;
  }

protected:
  std::int32_t native_ ;
  bool         managed_;
};
}