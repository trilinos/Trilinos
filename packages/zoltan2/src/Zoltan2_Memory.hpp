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


namespace Z2{

#ifdef Z2_OMIT_ALL_ERROR_CHECKING

#define Z2_SYNC_MEMORY_ALLOC(comm, env, datatype, ptrname, nbytes){ \
  datatype *ptrname = NULL; \
  if (nbytes) \
    ptrname = new datatype [nbytes]; \
}

#define Z2_ASYNC_MEMORY_ALLOC(comm, env, datatype, ptrname, nbytes){ \
  datatype *ptrname = NULL; \
  if (nbytes) \
    ptrname = new datatype [nbytes]; \
}

#else

/*! Allocate memory followed by a global check of success.

    All throw an error if any failed.
 */

#define Z2_SYNC_MEMORY_ALLOC(comm, env, datatype, ptrname, num) {\
  ptrname = NULL; \
  int fail = 0, gfail=0; \
  if ((num) > 0){ \
    ptrname = new datatype [num]; \
    if (!ptrname) fail = 1;  \
  } \
  Teuchos::reduceAll<int, int>(comm, Teuchos::REDUCE_MAX, 1, &fail, &gfail); \
  if (gfail > 0) { \
    if (fail > 0) \
      *(env)._errorOStream << (comm).getRank() << ": " << __FILE__ << ", " << __LINE__ << ", size " << num << std::endl; \
    throw std::bad_alloc(); \
  } \
}

/*! Allocate memory, throw error if it fails.
 */

#define Z2_ASYNC_MEMORY_ALLOC(comm, env, datatype, ptrname, num) { \
  ptrname = NULL; \
  if ((num) > 0) {\
    ptrname = new datatype [num]; \
    if (!ptrname) { \
      *(env)._errorOStream << (comm).getRank() << ": " << __FILE__ << ", " << __LINE__ << ", size " << num << std::endl; \
      throw std::bad_alloc(); \
    } \
  } \
}

#endif

} //namespace Z2

#endif


