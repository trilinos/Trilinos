// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

#ifndef _ZOLTAN2_EXECEPTIONS_HPP_
#define _ZOLTAN2_EXECEPTIONS_HPP_

/*! \file Zoltan2_Exceptions.hpp

  Exception handling macros.  We throw 3 types of error:

   \li \c  std::runtime_error   for an error in input
   \li \c  std::bad_alloc       for failure to allocate memory
   \li \c  std::logic_error     for an apparent bug in the code

  The GLOBAL macros are for assertions that all processes in
  a communicator call.  The throw an error if any of the
  processes finds that the assertion fails.

  The LOCAL macros are local.  If the assertion fails, the
  process throws an error.
  
*/

#include <stdexcept>
#include <iostream>
#include <string>
#include <Zoltan2.hpp>

/*!  We should always check basic assertions.
*/
#define Z2_BASIC_ASSERTION      0

/*!  Extra, more expensive level of checking.
 *
 * A parameter will state whether "extra" checking should be
 * done.  An example of extra checking is checking that an
 * input graph is value.  
 */
#define Z2_COMPLEX_ASSERTION    1

/*!  Even more extensive checking.
 *
 * This is extra checking we would do when debugging
 * a problem.
 *
 */
#define Z2_DEBUG_MODE_ASSERTION  2

#define Z2_LOCAL_INPUT_ASSERTION( \
Zoltan_Parameters params, s, assertion, level) \
{ \
  int check_level = params.check_level; \
  if (level <= check_level) {
    if (!(assertion)){ \
      ostringstream oss; \
      oss << ___FILE___ << ":" << __LINE__; \
      if (s.size() > 0) oss << ": " <<s; \
      oss << std::endl; \
      throw(std::runtime_error(oss.str()); \
    } \
  } \
}
#define Z2_LOCAL_MEMORY_ASSERTION( \
Zoltan_Parameters params, allocSize, assertion) \
{ \
  if (!(assertion)){ \
    ostringstream oss; \
    oss << ___FILE___ << ":" << __LINE__ << " " << (allocSize) << std::endl; \
    throw(std::bad_alloc(oss.str()); \
  } \
}

#define Z2_LOCAL_BUG_ASSERTION( \
Zoltan_Parameters params, s, assertion, level) \
{ \
  int check_level = params.check_level; \
  if (level <= check_level) {
    if (!(assertion)){ \
      ostringstream oss; \
      oss << ___FILE___ << ":" << __LINE__; \
      if (s.size() > 0) oss << ": " <<s; \
      oss << std::endl; \
      throw(std::logic_error(oss.str()); \
    } \
  } \
}

#define Z2_GLOBAL_INPUT_ASSERTION( \
MPI_Communicator comm, Zoltan_Parameters params, s, assertion, level) \
{ \
  int check_level = params.check_level; \
  if (level <= check_level) {  \
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

#define Z2_GLOBAL_MEMORY_ASSERTION( \
MPI_Communicator comm, Zoltan_Parameters params, allocSize, assertion) \
{ \
  int pass = 1, gpass=0;  \
  if (!(assertion)) pass = 0;  \
  MPI_Allreduce(&pass, &gpass, 1, MPI_INT, MPI_MAX, comm);  \
  if (gpass > 0){  \
    ostringstream oss; \
    oss << ___FILE___ << ":" << __LINE__ << " " << (allocSize) << std::endl;\
    throw(std::bad_alloc(oss.str()); \
  } \
}

#define Z2_GLOBAL_BUG_ASSERTION( \
MPI_Communicator comm, Zoltan_Parameters params, s, assertion, level) \
{ \
  int check_level = params.check_level; \
  if (level <= check_level) {  \
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

/*! Throw an error returned from outside the Zoltan2 library.
 */
#define Z2_THROW_OUTSIDE_ERROR(Zoltan_Parameters params, std::exception e) \
{ \
  ostream &os = params.debug_stream; \
  os << ___FILE___ << ":" << __LINE__ << ": " << e.what() << std::endl; \
  throw(e); \
}
   
#endif
