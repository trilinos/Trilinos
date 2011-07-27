// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

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

/*!  We should always check basic assertions.
*/
#define Z2_BASIC_ASSERTION      0

/*!  Extra, more expensive level of checking.
 *
 * A parameter will state whether "extra" checking should be
 * done.  An example of extra checking is checking that an
 * input graph is valid.  
 */
#define Z2_COMPLEX_ASSERTION    1

/*!  Even more extensive checking.
 *
 * This is extra checking we would do when debugging
 * a problem.
 *
 */
#define Z2_DEBUG_MODE_ASSERTION  2

#define Z2_MAX_CHECK_LEVEL Z2_DEBUG_MODE_ASSERTION  

#ifdef ZOLTAN2_OMIT_ALL_ERROR_CHECKING

#define Z2_LOCAL_INPUT_ASSERTION(env, s, assertion, level) {}
#define Z2_LOCAL_BUG_ASSERTION( env, s, assertion, level) {}
#define Z2_LOCAL_MEMORY_ASSERTION( env, requestSize, assertion) {}
#define Z2_GLOBAL_INPUT_ASSERTION( comm, env, s, assertion, level) {}
#define Z2_GLOBAL_BUG_ASSERTION( comm, env, s, assertion, level) {}
#define Z2_GLOBAL_MEMORY_ASSERTION( comm, env, requestSize, assertion) {}

#else

#define Z2_LOCAL_INPUT_ASSERTION(env, s, assertion, level) { \
  if (level <= env.errorCheckLevel) {
    if (!(assertion)){ \
      ostringstream oss; \
      oss << ___FILE___ << ":" << __LINE__; \
      if (s.size() > 0) oss << ": " <<s; \
      oss << std::endl; \
      throw(std::runtime_error(oss.str()); \
    } \
  } \
}

#define Z2_LOCAL_BUG_ASSERTION( env, s, assertion, level) { \
  if (level <= env.errorCheckLevel) {
    if (!(assertion)){ \
      ostringstream oss; \
      oss << ___FILE___ << ":" << __LINE__; \
      if (s.size() > 0) oss << ": " <<s; \
      oss << std::endl; \
      throw(std::logic_error(oss.str()); \
    } \
  } \
}

#define Z2_LOCAL_MEMORY_ASSERTION( env, requestSize, assertion) {
  if (!(assertion)){ \
    ostringstream oss; \
    oss << ___FILE___ << ":" << __LINE__ << " size " << requestSize; \
    oss << std::endl; \
    throw(std::bad_alloc(oss.str()); \
  } \
}

#define Z2_GLOBAL_INPUT_ASSERTION( comm, env, s, assertion, level) { \
  if (level <= env.errorCheckLevel) {  \
    int pass = 1, gpass=0;  \
    if (!(assertion)) pass = 0;  \
    MPI_Allreduce(&pass, &gpass, 1, MPI_INT, MPI_MAX, comm);  \
    if (gpass > 0){  \
      ostringstream oss; \
      oss << ___FILE___ << ":" << __LINE__; \
      if (s.size() > 0) oss << ": " <<s; \
      oss << std::endl; \
      throw(std::runtime_error(oss.str()); \
    } \
  } \
}

#define Z2_GLOBAL_BUG_ASSERTION( comm, env, s, assertion, level) \
{ \
  int env.errorCheckLevel = env.check_level; \
  if (level <= env.errorCheckLevel) {  \
    int pass = 1, gpass=0;  \
    if (!(assertion)) pass = 0;  \
    MPI_Allreduce(&pass, &gpass, 1, MPI_INT, MPI_MAX, comm);  \
    if (gpass > 0){  \
      ostringstream oss; \
      oss << ___FILE___ << ":" << __LINE__; \
      if (s.size() > 0) oss << ": " <<s; \
      oss << std::endl; \
      throw(std::logic_error(oss.str()); \
    } \
  } \
}

#define Z2_GLOBAL_MEMORY_ASSERTION( comm, env, requestSize, assertion) \
{\
  int pass = 1, gpass=0;  \
  if (!(assertion)) pass = 0;  \
  MPI_Allreduce(&pass, &gpass, 1, MPI_INT, MPI_MAX, comm);  \
  ostringstream oss; \
  if (pass > 0){  \
    oss << ___FILE___ << ":" << __LINE__ << "size " << requestSize; \
    oss << std::endl; \
  } \
  if (gpass > 0)  \
    throw(std::bad_alloc(oss.str()); \
}

#endif

/*! Throw an error returned from outside the Zoltan2 library.
 */
#define Z2_THROW_OUTSIDE_ERROR(env, e) \
{ \
  ostream &os = env.debug_stream; \
  os << ___FILE___ << ":" << __LINE__; \
  if (e.what.size() > 0) \
    os << ": " << e.what();  \
  os << std::endl; \
  throw(e); \
}

/*! Throw an error returned from another Zoltan2 method.
 */
#define Z2_THROW_ZOLTAN2_ERROR(env, e) { throw(e); }
   
#endif
