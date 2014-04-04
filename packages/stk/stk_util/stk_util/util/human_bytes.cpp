#include <stk_util/util/human_bytes.hpp>

#include <string>
#include <sstream>
#include <unistd.h>
#include <iomanip>

namespace stk
{
  std::string human_bytes(size_t arg_bytes)
  {
    double bytes = arg_bytes;
    const double K = 1024;
    const double M = K*1024;
    const double G = M*1024;

    std::ostringstream out;
    if (bytes < K) {
      out << std::setprecision(4) << bytes << " B";
    } else if (bytes < M) {
      bytes /= K;
      out << std::setprecision(4) << bytes << " K";
    } else if (bytes < G) {
      bytes /= M;
      out << std::setprecision(4) << bytes << " M";
    } else {
      bytes /= G;
      out << std::setprecision(4) << bytes << " G";
    }
    return out.str();
  }
}
