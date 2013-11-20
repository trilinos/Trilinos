#ifndef STK_IO_UTILS_h
#define STK_IO_UTILS_h

#include <string>

namespace stk {
  namespace io {
    class Utils {
    public:

      static void check_shutdown_request();
      static unsigned int check_control_request(double time, int step);

      /*! See if filename contains "%P" which is replaced by the number of processors...
       * Assumes that %P only occurs once...
       * filename is changed.
       */
      static void filename_substitution(std::string &filename);
    };
  }
}
#endif
