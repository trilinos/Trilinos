#ifndef KRINO_KRINO_MESH_UTILS_AKRI_ALLREDUCE_HPP_
#define KRINO_KRINO_MESH_UTILS_AKRI_ALLREDUCE_HPP_

#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/parallel/Parallel.hpp>

namespace krino {

template <typename... Args> std::array<typename std::common_type<Args...>::type, sizeof...(Args)> make_array(Args&&... args)
{
  return {{ args... }};
}

template<typename T, size_t N>
std::array<T,N> all_reduce_sum_array(stk::ParallelMachine comm, const std::array<T,N> & local)
{
  std::array<T,N>  global;
  stk::all_reduce_sum(comm, local.data(), global.data(), local.size());
  return global;
}

template<typename T, size_t N>
std::array<T,N> all_reduce_min_array(stk::ParallelMachine comm, const std::array<T,N> & local)
{
  std::array<T,N>  global;
  stk::all_reduce_min(comm, local.data(), global.data(), local.size());
  return global;
}

template<typename T, size_t N>
std::array<T,N> all_reduce_max_array(stk::ParallelMachine comm, const std::array<T,N> & local)
{
  std::array<T,N>  global;
  stk::all_reduce_max(comm, local.data(), global.data(), local.size());
  return global;
}

template<typename... Args>
void all_reduce_sum(stk::ParallelMachine comm, Args&&... args)
{
  const auto global = all_reduce_sum_array(comm, make_array(args...));
  const auto * data = global.data();
  ((std::forward<Args>(args) = *(data++)), ...);
}

template<typename... Args>
void all_reduce_min(stk::ParallelMachine comm, Args&&... args)
{
  const auto global = all_reduce_min_array(comm, make_array(args...));
  const auto * data = global.data();
  ((std::forward<Args>(args) = *(data++)), ...);
}

template<typename... Args>
void all_reduce_max(stk::ParallelMachine comm, Args&&... args)
{
  const auto global = all_reduce_max_array(comm, make_array(args...));
  const auto * data = global.data();
  ((std::forward<Args>(args) = *(data++)), ...);
}

}



#endif /* KRINO_KRINO_MESH_UTILS_AKRI_ALLREDUCE_HPP_ */
