#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <functional>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#ifdef ASTRAY_USE_MPI
#include <mpi.h>
#endif

namespace ast
{
template <typename type>
struct record
{
  type        mean              ()                             
  {
    return std::accumulate(values.begin(), values.end(), type(0)) / values.size();
  }
  type        variance          ()                            
  {
    auto m = mean();
    std::vector<type> differences(values.size());
    std::transform(values.begin(), values.end(), differences.begin(), [m] (const type& value) { return value - m; });
    return std::inner_product(differences.begin(), differences.end(), differences.begin(), type(0)) / values.size();
  }
  type        standard_deviation()                            
  {
    return std::sqrt(variance());
  }

  std::string to_string         ()                            
  {
    std::ostringstream stream;
    stream.precision(std::numeric_limits<type>::max_digits10);
    stream << name << ",";
    for (auto& value : values)
      stream << value << ",";
    stream << mean() << "," << variance() << "," << standard_deviation();
    return stream.str();
  }
  void        to_csv            (const std::string& filepath) 
  {
    std::ofstream stream(filepath);
    stream << "name,";
    for (auto i = 0; i < values.size(); ++i)
      stream << "run_" << i << ",";
    stream << "mean,variance,standard deviation\n";
    stream << to_string();
  }

  std::string       name  ;
  std::vector<type> values;
};

template <typename type>
struct session
{
  session           ()                     = default;
  session           (const session&  that) = default;
  session           (      session&& temp) = default;
  virtual ~session  ()                     = default;
  session& operator=(const session&  that) = default;
  session& operator=(      session&& temp) = default;

  virtual std::string to_string()
  {
    std::ostringstream stream;
    for (auto& record : records)
      stream << record.to_string() << "\n";
    return stream.str();
  }
  virtual void        to_csv   (const std::string& filepath)
  {
    std::ofstream stream(filepath);
    stream << "name,";
    for (auto i = 0; i < records[0].values.size(); ++i)
      stream << "run_" << i << ",";
    stream << "mean,variance,standard deviation\n";
    stream << to_string();
  }

  std::vector<record<type>> records;
};

#ifdef ASTRAY_USE_MPI
template <typename type>
class  mpi_session : public session<type>
{
public:
  explicit mpi_session  (const MPI_Comm communicator, const std::int32_t master_rank) : communicator_(communicator), master_rank_(master_rank)
  {
    MPI_Comm_rank(communicator_, &rank_);
    MPI_Comm_size(communicator_, &size_);
  }
  mpi_session           (const mpi_session&  that) = default;
  mpi_session           (      mpi_session&& temp) = default;
 ~mpi_session           () override                = default;
  mpi_session& operator=(const mpi_session&  that) = default;
  mpi_session& operator=(      mpi_session&& temp) = default;
  
  void        gather   ()
  {
    std::ostringstream stream;
    for (auto& record : session<type>::records)
      stream << rank_ << "," << record.to_string() << "\n";
    const std::string local_string = stream.str();
    const auto        local_size   = static_cast<std::int32_t>(local_string.size());
    
    std::vector<std::int32_t> sizes        (size_);
    std::vector<std::int32_t> displacements(size_);
    std::int32_t              counter = 0;
    MPI_Gather (&local_size, 1, MPI_INT, sizes.data(), 1, MPI_INT, master_rank_, communicator_);
    for (auto i = 0; i < size_; ++i)
      displacements[i] = counter, counter += sizes[i];
    gathered_.resize(counter);
    MPI_Gatherv(local_string.data(), local_size, MPI_CHAR, gathered_.data(), sizes.data(), displacements.data(), MPI_CHAR, master_rank_, communicator_);
  }
  std::string to_string()                            override
  {
    return rank_ == master_rank_ ? gathered_ : session<type>::to_string();
  }
  void        to_csv   (const std::string& filepath) override
  {
    if (rank_ != master_rank_)
      return;

    std::ofstream stream(filepath);
    stream << "rank,name,";
    for (auto i = 0; i < session<type>::records[0].values.size(); ++i)
      stream << "run_" << i << ",";
    stream << "mean,variance,standard deviation\n";
    stream << to_string();
  }
  
protected:
  MPI_Comm     communicator_;
  std::int32_t master_rank_ ;
  std::int32_t rank_        ;
  std::int32_t size_        ;
  std::string  gathered_    ;
};
#endif

template <typename type, typename period, typename clock>
class  session_recorder
{
public:
  explicit session_recorder  (const std::size_t index, const std::size_t iterations, session<type>& session) 
  : index_(index), iterations_(iterations), session_(session)
  {

  }
  session_recorder           (const session_recorder&  that) = delete ;
  session_recorder           (      session_recorder&& temp) = default;
  virtual ~session_recorder  ()                              = default;
  session_recorder& operator=(const session_recorder&  that) = delete ;
  session_recorder& operator=(      session_recorder&& temp) = default;
  
  void record(const std::string& name, const std::function<void()>& function)
  {
    const auto start = clock::now();
    function();
    const auto end   = clock::now();

    auto record = std::find_if(session_.records.begin(), session_.records.end(),
      [&name] (const ast::record<type>& iteratee) { return iteratee.name == name; });
    if (record == session_.records.end())
    {
      session_.records.push_back({name, {std::vector<type>(iterations_)}});
      record = std::prev(session_.records.end());
    }
    record->values[index_] = std::chrono::duration<type, period>(end - start).count();
  }

protected:
  const std::size_t index_     ;
  const std::size_t iterations_;
  session<type>&    session_   ;
};

template<typename type, typename period, typename clock>
record<type>      benchmark    (const std::function<void()>&                                       function, const std::size_t iterations = 1)
{
  record<type> record {"benchmark", std::vector<type>(iterations)};
  for (auto i = 0; i < iterations; ++i)
  {
    const auto start = clock::now();
    function();
    const auto end   = clock::now();
    record.values[i] = std::chrono::duration<type, period>(end - start).count();
  }
  return record;
}
template<typename type, typename period, typename clock>
session<type>     benchmark    (const std::function<void(session_recorder<type, period, clock>&)>& function, const std::size_t iterations = 1)
{
  session<type> session;
  for(auto i = 0; i < iterations; ++i)
  {
    session_recorder<type, period, clock> recorder(i, iterations, session);
    function(recorder);
  }
  return session;
}

#ifdef ASTRAY_USE_MPI
template<typename type, typename period, typename clock>
mpi_session<type> benchmark_mpi(const std::function<void(session_recorder<type, period, clock>&)>& function, const std::size_t iterations = 1, 
  const MPI_Comm communicator = MPI_COMM_WORLD, const std::int32_t master_rank = 0)
{
  mpi_session<type> session(communicator, master_rank);
  for (auto i = 0; i < iterations; ++i)
  {
    session_recorder<type, period, clock> recorder(i, iterations, session);
    function(recorder);
  }
  session.gather();
  return session;
}
#endif
}