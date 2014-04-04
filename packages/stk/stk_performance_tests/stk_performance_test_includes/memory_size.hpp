#ifndef samba_performance_memory_size_hpp
#define samba_performance_memory_size_hpp

#include <iostream>
#include <stk_util/util/human_bytes.hpp>
#include <stk_util/environment/memory_util.hpp>

namespace samba { namespace detail {

inline
void report_memory_size(std::ostream & out, const std::string & tag = std::string(""))
{
  size_t now = 0;
  size_t hwm = 0;
  stk::get_memory_usage(now, hwm);

  out << tag << " Current Memory Used:  " << now << " (" << stk::human_bytes(now) << ")" << std::endl;
  out << tag << "     Max Memory Used:  " << hwm << " (" << stk::human_bytes(hwm) << ")" << std::endl;
}

} } // namespace samba::detail


#endif // samba_performance_memory_size_hpp
