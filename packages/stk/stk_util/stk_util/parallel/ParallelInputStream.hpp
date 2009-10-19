#ifndef stk_util_parallel_ParallelInputStream_hpp
#define stk_util_parallel_ParallelInputStream_hpp

#include <stk_util/parallel/Parallel.hpp>
#include <istream>

namespace stk {

class ParallelInputStream : public std::istream {
public:
  ParallelInputStream( ParallelMachine comm , const char * const file_name );
  virtual ~ParallelInputStream();
};

}

//----------------------------------------------------------------------

#endif

