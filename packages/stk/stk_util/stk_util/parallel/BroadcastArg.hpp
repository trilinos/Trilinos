#ifndef STK_UTIL_PARALLEL_BROADCASTARG_HPP
#define STK_UTIL_PARALLEL_BROADCASTARG_HPP

#include <stk_util/parallel/Parallel.hpp>

namespace stk {

/**
 * @brief Class <b>BroadcastArg</b> creates a copy of argc and argv after broadcasting
 * them from processor 0.
 *
 */
struct BroadcastArg 
{
  /**
   * Creates a new <b>BroadcastArg</b> instance.
   *
   * @param parallel_machine	        a <b>ParallelMachine</b> value that species the
   *                                    parallel machine from performing the broadcast.
   *
   * @param argc			an <b>int</b> value of argc from main
   *
   * @param argv			a <b>char</b> pointer pointer of argv from main
   *
   */
  BroadcastArg(ParallelMachine parallel_machine, int argc, char **argv);

  /**
   * Destroys a <b>BroadcastArg</b> instance.
   *
   */
  ~BroadcastArg();

  int           m_argc;                         ///< The broadcasted argc
  char **       m_argv;                         ///< The broadcasted argv

private:
  BroadcastArg(const BroadcastArg &argv);
  BroadcastArg &operator=(const BroadcastArg &argv);
};

} // namespace stk

#endif // STK_UTIL_PARALLEL_BROADCASTARG_HPP
