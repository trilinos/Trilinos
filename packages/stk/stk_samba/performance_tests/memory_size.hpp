#ifndef samba_performance_memory_size_hpp
#define samba_performance_memory_size_hpp

#include <iostream>
#include <sstream>

#include <stk_util/util/MallocUsed.h>

namespace samba { namespace detail {

inline
std::string human_bytes(double bytes)
{
  const size_t K = 1024;
  const size_t M = K*1024;
  const size_t G = M*1024;

  std::ostringstream out;
  if (bytes < K) {
    out << bytes << " B";
  } else if (bytes < M) {
    bytes /= K;
    out << bytes << " K";
  } else if (bytes < G) {
    bytes /= M;
    out << bytes << " M";
  } else {
    bytes /= G;
    out << bytes << " G";
  }
  return out.str();
}

inline
void report_memory_size(std::ostream & out, const std::string & tag = std::string(""))
{
  out << tag << "          Current Memory Used:  " << malloc_used() << " (" << human_bytes(malloc_used()) << ")" << std::endl;
  out << tag << "     Current Memory Footprint:  " << malloc_footprint() << " (" << human_bytes(malloc_footprint()) << ")" << std::endl;
  out << tag << " Current Max Memory Footprint:  " << malloc_max_footprint() << " (" << human_bytes(malloc_max_footprint()) << ")" << std::endl;
}

} } // namespace samba::detail


#endif // samba_performance_memory_size_hpp
