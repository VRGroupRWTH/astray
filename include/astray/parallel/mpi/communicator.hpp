#pragma once

#include <cstdint>

#ifdef ASTRAY_USE_MPI
#include <mpi.h>
#endif

#include <astray/parallel/mpi/data_type.hpp>

namespace ast::mpi
{
class communicator
{
public:
#ifdef ASTRAY_USE_MPI
  explicit communicator  (const std::int32_t native = MPI_COMM_WORLD)
#else
  explicit communicator  (const std::int32_t native = 0)
#endif
  : native_(native)
  {
  
  }
  communicator           (const communicator&  that) = default;
  communicator           (      communicator&& temp) = default;
 ~communicator           ()                          = default;
  communicator& operator=(const communicator&  that) = default;
  communicator& operator=(      communicator&& temp) = default;
  
  template <typename type>
  void                   bcast  (type* data, std::int32_t size, const data_type& data_type, const std::int32_t root = 0)
  {
#ifdef ASTRAY_USE_MPI
    MPI_Bcast(data, size, data_type.native(), root, native_);
#endif
  }
  
  template <typename type>
  void                   gatherv(
    type* sent    , std::int32_t  sent_size     ,                              const data_type& sent_data_type    ,
    type* received, std::int32_t* received_sizes, std::int32_t* displacements, const data_type& received_data_type, const std::int32_t root = 0)
  {
#ifdef ASTRAY_USE_MPI
    MPI_Gatherv(sent, sent_size, sent_data_type.native(), received, received_sizes, displacements, received_data_type.native(), root, native_);
#endif
  }

  std::int32_t           rank   () const
  {
    std::int32_t result(0);
#ifdef ASTRAY_USE_MPI
    MPI_Comm_rank(native_, &result);
#endif
    return result;
  }
  std::int32_t           size   () const
  {
    std::int32_t result(1);
#ifdef ASTRAY_USE_MPI
    MPI_Comm_size(native_, &result);
#endif
    return result;
  }

  constexpr std::int32_t native () const
  {
    return native_;
  }

protected:
  std::int32_t native_;
};
}