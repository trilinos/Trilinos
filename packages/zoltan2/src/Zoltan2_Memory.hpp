// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

#ifndef _ZOLTAN2_MEMORY_HPP_
#define _ZOLTAN2_MEMORY_HPP_

/*! \file Zoltan2_Memory.hpp

  Memory allocation macros.

  If we wanted to do clever memory management, we would do it here.  This
  would probably require the Zoltan2::Environment object, which is why
  it is passed.
*/

#include <Teuchos_CommHelpers.hpp>
#include <Zoltan2_config.h>


namespace Zoltan2{

#ifdef HAVE_MALLINFO
int getAllocatedMemory();
int mallocCount(std::string label);
int getMallocCount(std::string label);
void printMallocCount();
void eraseMallocCount();
#endif

#ifdef Z2_OMIT_ALL_ERROR_CHECKING

#define Z2_SYNC_MEMORY_ALLOC(comm, env, datatype, ptrname, num){ \
  datatype *ptrname = NULL; \
  if ((num) > 1) \
    ptrname = new datatype [num]; \
  else if ((num) > 0) \
    ptrname = new datatype; \
}

#define Z2_ASYNC_MEMORY_ALLOC(comm, env, datatype, ptrname, num){ \
  datatype *ptrname = NULL; \
  if ((num) > 1) \
    ptrname = new datatype [num]; \
  else if ((num) > 0) \
    ptrname = new datatype; \
}

#else

/*! Allocate memory followed by a global check of success.

    All throw an error if any failed.
 */

#define Z2_SYNC_MEMORY_ALLOC(comm, env, datatype, ptrname, num) {\
  ptrname = NULL; \
  int fail = 0, gfail=0; \
  if ((num) > 1){ \
    ptrname = new datatype [num]; \
    if (!ptrname) fail = 1;  \
  } \
  else if ((num) == 1) { \
    ptrname = new datatype; \
    if (!ptrname) fail = 1;  \
  } \
  Teuchos::reduceAll<int, int>(comm, Teuchos::REDUCE_MAX, 1, &fail, &gfail); \
  if (gfail > 0) { \
    if (fail > 0){ \
      std::ostringstream _msg;  \
      _msg << __FILE__ << ", " << __LINE__ << ", size " << num << std::endl; \
      (env).dbg_->error(_msg.str()); \
    } \
    throw std::bad_alloc(); \
  } \
}

/*! Allocate memory, throw error if it fails.
 */

#define Z2_ASYNC_MEMORY_ALLOC(comm, env, datatype, ptrname, num) { \
  ptrname = NULL; \
  if ((num) > 1) {\
    ptrname = new datatype [num]; \
  } \
  else if ((num) == 1){ \
    ptrname = new datatype; \
  } \
  if ((num) && !ptrname) { \
    std::ostringstream _msg;  \
    _msg << __FILE__ << ", " << __LINE__ << ", size " << num << std::endl; \
    (env).dbg_->error(_msg.str()); \
    throw std::bad_alloc(); \
  } \
}

#endif

} //namespace Zoltan2

#endif


