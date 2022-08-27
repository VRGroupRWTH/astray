#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>

#include <astray/math/indexing.hpp>
#include <astray/math/prime_factorization.hpp>

namespace ast
{
// Partitions an N-dimensional domain of domain_size to a hyper-rectangular grid of grid_size, 
// which consists of blocks of block_size based on the communicator_size.
// Also computes a rank_multi_index which is the N-dimensional index of the rank, and a 
// rank_offset which equals block_size * rank_multi_index, based on the communicator_rank.
template <std::size_t dimensions, typename size_type = std::size_t, typename multi_size_type = std::array<size_type, dimensions>, bool fortran_order = false>
class partitioner
{
public:
  explicit partitioner(
    const size_type        communicator_rank, 
    const size_type        communicator_size, 
    const multi_size_type& domain_size      )
  : communicator_rank_(communicator_rank)
  , communicator_size_(communicator_size)
  , domain_size_      (domain_size      )
  {
    partitioner<dimensions, size_type, multi_size_type, fortran_order>::update();
  }
  partitioner           (const partitioner&  that) = default;
  partitioner           (      partitioner&& temp) = default;
  virtual ~partitioner  ()                         = default;
  partitioner& operator=(const partitioner&  that) = default;
  partitioner& operator=(      partitioner&& temp) = default;

  void                   set_communicator_rank(const size_type        communicator_rank)
  {
    communicator_rank_ = communicator_rank;
    update();
  }
  void                   set_communicator_size(const size_type        communicator_size)
  {
    communicator_size_ = communicator_size;
    update();
  }
  void                   set_domain_size      (const multi_size_type& domain_size      )
  {
    domain_size_ = domain_size;
    update();
  }
  
  size_type              communicator_rank    () const
  {
    return communicator_rank_;
  }
  size_type              communicator_size    () const
  {
    return communicator_size_;
  }
  const multi_size_type& domain_size          () const
  {
    return domain_size_;
  }

  const multi_size_type& grid_size            () const
  {
    return grid_size_;
  }
  const multi_size_type& block_size           () const
  {
    return block_size_;
  }
  const multi_size_type& rank_multi_index     () const
  {
    return rank_multi_index_;
  }
  const multi_size_type& rank_offset          () const
  {
    return rank_offset_;
  }

protected:
  virtual void           update               ()
  {
    auto prime_factors = prime_factorize(communicator_size_);
    auto current_size  = domain_size_;
    grid_size_.fill(1);
    while (!prime_factors.empty())
    {
      auto dimension = std::distance(current_size.begin(), std::max_element(current_size.begin(), current_size.end()));
      current_size[dimension] /= prime_factors.back();
      grid_size_  [dimension] *= prime_factors.back();
      prime_factors.pop_back();
    }
    std::transform(domain_size_.begin(), domain_size_.end(), grid_size_       .begin(), block_size_ .begin(), std::divides   <>());
    rank_multi_index_ = unravel_index<multi_size_type, fortran_order>(communicator_rank_, grid_size_);
    std::transform(block_size_ .begin(), block_size_ .end(), rank_multi_index_.begin(), rank_offset_.begin(), std::multiplies<>());
  }
  
  size_type       communicator_rank_;
  size_type       communicator_size_;
  multi_size_type domain_size_      ;

  multi_size_type grid_size_        ;
  multi_size_type block_size_       ;
  multi_size_type rank_multi_index_ ;
  multi_size_type rank_offset_      ;
};
}