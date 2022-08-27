#pragma once

#ifdef ASTRAY_USE_MPI
#include <mpi.h>
#endif

namespace ast::mpi
{
class environment
{
public:
  environment           ()
  {
#ifdef ASTRAY_USE_MPI
    auto initialized = 0;
    MPI_Initialized(&initialized);
    if (!initialized)
    {
      MPI_Init(nullptr, nullptr);
      finalize_ = true;
    }
#endif
  }
  environment           (const environment&  that) = delete;
  environment           (      environment&& temp) = delete;
  virtual ~environment  ()
  {
#ifdef ASTRAY_USE_MPI
    auto finalized = 0;
    MPI_Finalized(&finalized);
    if (finalize_ && !finalized)
      MPI_Finalize();
#endif
  }
  environment& operator=(const environment&  that) = delete;
  environment& operator=(      environment&& temp) = delete;

protected:
  bool finalize_ = false;
};
}