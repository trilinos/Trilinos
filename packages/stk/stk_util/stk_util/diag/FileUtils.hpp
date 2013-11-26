#ifndef STK_util_FILEUTILS_h
#define STK_util_FILEUTILS_h

#include <string>

namespace stk {
  namespace util {
    void check_shutdown_request();
    unsigned int check_control_request(double time, int step);

      /*! See if filename contains "%P" which is replaced by the number of processors...
       * Assumes that %P only occurs once...
       * filename is changed.
       */
    void filename_substitution(std::string &filename);
  }
}
#endif
