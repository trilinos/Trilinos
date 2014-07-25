/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_util_parallel_ParallelInputStream_hpp
#define stk_util_parallel_ParallelInputStream_hpp

#include <stk_util/parallel/Parallel.hpp>
#include <istream>

namespace stk_classic {

class ParallelInputStream : public std::istream {
public:
  ParallelInputStream( ParallelMachine comm , const char * const file_name );
  virtual ~ParallelInputStream();
};

}

//----------------------------------------------------------------------

#endif

