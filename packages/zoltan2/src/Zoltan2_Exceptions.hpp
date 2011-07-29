// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_Exceptions.hpp

    \brief Exception handling macros
*/


#ifndef _ZOLTAN2_EXCEPTIONS_HPP_
#define _ZOLTAN2_EXCEPTIONS_HPP_

/*! \file Zoltan2_Exceptions.hpp

  Exception handling macros.  We throw 3 types of error:

   \li \c  std::runtime_error   for an error in input
   \li \c  std::bad_alloc       for failure to allocate memory
   \li \c  std::logic_error     for an apparent bug in the code

  The GLOBAL macros are for assertions that all processes in
  a communicator call.  They throw an error if any of the
  processes finds that the assertion fails.

  The LOCAL macros are local.  If the assertion fails, the
  process throws an error.

  The bad_alloc exceptions are thrown in Zoltan2_Memory.hpp.
*/

#include <stdexcept>
#include <iostream>
#include <string>
#include <Zoltan2_Parameters.hpp>

#ifdef Z2_OMIT_ALL_ERROR_CHECKING

#define Z2_LOCAL_INPUT_ASSERTION(comm, env, s, assertion, level) {}
#define Z2_LOCAL_BUG_ASSERTION(comm,  env, s, assertion, level) {}
#define Z2_LOCAL_MEMORY_ASSERTION(comm,  env, requestSize, assertion) {}
#define Z2_GLOBAL_INPUT_ASSERTION( comm, env, s, assertion, level) {}
#define Z2_GLOBAL_BUG_ASSERTION( comm, env, s, assertion, level) {}
#define Z2_GLOBAL_MEMORY_ASSERTION( comm, env, requestSize, assertion) {}

#else

#define Z2_LOCAL_INPUT_ASSERTION(comm, env, s, assertion, level) { \
  if (level <= (env)._errorCheckLevel) { \
    if (!(assertion)){ \
      std::ostringstream oss; \
      oss << (comm).getRank() << ": " << __FILE__ << ", " << __LINE__ << ", " << s << std::endl; \
      throw std::runtime_error(oss.str()); \
    } \
  } \
}

#define Z2_LOCAL_BUG_ASSERTION( comm, env, s, assertion, level) { \
  if (level <= (env)._errorCheckLevel) { \
    if (!(assertion)){ \
      std::ostringstream oss; \
      oss << (comm).getRank() << ": " << __FILE__ << ", " << __LINE__ << ", " << s << std::endl; \
      throw std::logic_error(oss.str()); \
    } \
  } \
}

/*! We always check for success of memory allocation, regardless of ERROR_CHECK_LEVEL.
 */

#define Z2_LOCAL_MEMORY_ASSERTION( comm, env, requestSize, assertion) { \
  if (!(assertion)){ \
    *(env)._errorOStream << (comm).getRank() << ": " << __FILE__ << ", " << __LINE__ << ", size " << requestSize << std::endl; \
    throw std::bad_alloc(); \
  } \
}

#define Z2_GLOBAL_INPUT_ASSERTION( comm, env, s, assertion, level) { \
  if (level <= (env)._errorCheckLevel) {  \
    int fail = 0, gfail=0;  \
    if (!(assertion)) fail = 1;  \
    Teuchos::reduceAll<int, int>(comm, Teuchos::REDUCE_MAX, 1, &fail, &gfail); \
    if (gfail > 0){  \
      std::ostringstream oss; \
      if (fail > 0) \
        oss << (comm).getRank() << ": " << __FILE__ << ", " << __LINE__ << ", " << s << std::endl; \
      throw std::runtime_error(oss.str()); \
    } \
  } \
}

#define Z2_GLOBAL_BUG_ASSERTION( comm, env, s, assertion, level) { \
  if (level <= (env)._errorCheckLevel) {  \
    int fail = 0, gfail=0;  \
    if (!(assertion)) fail = 1;  \
    Teuchos::reduceAll<int, int>(comm, Teuchos::REDUCE_MAX, 1, &fail, &gfail); \
    if (gfail > 0){  \
      std::ostringstream oss; \
      if (fail > 0) \
        oss << (comm).getRank() << ": " << __FILE__ << ", " << __LINE__ << ", " << s << std::endl; \
      throw std::logic_error(oss.str()); \
    } \
  } \
}

/*! We always check for success of memory allocation, regardless of ERROR_CHECK_LEVEL.
 */

#define Z2_GLOBAL_MEMORY_ASSERTION( comm, env, requestSize, assertion) {\
  int fail = 0, gfail=0;  \
  if (!(assertion)) fail = 1;  \
  Teuchos::reduceAll<int, int>(comm, Teuchos::REDUCE_MAX, 1, &fail, &gfail); \
  if (gfail > 0){  \
    if (fail > 0) \
      *(env)._errorOStream << (comm).getRank() << ": " << __FILE__ << ", " << __LINE__ << ", size " << requestSize << std::endl; \
    throw std::bad_alloc(); \
  } \
}

#endif

/*! Throw an error returned from outside the Zoltan2 library.
 */
#define Z2_THROW_OUTSIDE_ERROR(env, e) { \
  *(env)._errorOStream << __FILE__ << ":" << __LINE__ << " " << e.what() << std::endl; \
  throw e; \
}

/*! Throw an error returned from another Zoltan2 method.
 */
#define Z2_THROW_ZOLTAN2_ERROR(env, e) { throw e; }
   
#endif
