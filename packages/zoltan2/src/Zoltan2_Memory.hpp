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

  If we wanted to do clever memory management, we would do it here.
*/

namespace Z2{

#ifdef ZOLTAN2_OMIT_ALL_ERROR_CHECKING

#define Z2_SYNC_MEMORY_ALLOC(comm, env, datatype, ptrname, nbytes) \
  datatype *ptrname = NULL; \
  if (nbytes) \
    ptrname = new datatype [nbytes]; 

#define Z2_ASYNC_MEMORY_ALLOC(comm, env, datatype, ptrname, nbytes) \
  datatype *ptrname = NULL; \
  if (nbytes) \
    ptrname = new datatype [nbytes]; 

#else

/*! Allocate memory followed by a global check of success.

    All throw an error if any failed.
 */

#define Z2_SYNC_MEMORY_ALLOC(comm, env, datatype, ptrname, num) \
  datatype *ptrname = NULL; \
  if (num) \
    ptrname = new datatype [num]; \
  int fail = 0, gfail=0; \
  if (num>0 && !ptrname) fail = 1;  \
  MPI_Allreduce(&fail, &gfail, 1, MPI_INT, MPI_MAX, comm);  \
  if (gfail > 0) { \
    ostringstream oss; \
    if (fail > 0){  \
      oss << "(" << comm.getRank() << ") " << ___FILE___ << ":" << __LINE__ << " " << (num) << std::endl;\
    } \
    throw(std::bad_alloc(oss.str()); \
  } \

/*! Allocate memory, throw error if it fails.
 */

#define Z2_ASYNC_MEMORY_ALLOC(comm, env, datatype, ptrname, num) \
{ \
  datatype *ptrname = NULL; \
  if (num > 0) {\
    ptrname = new datatype [num]; \
    if (!ptrname) { \
      ostringstream oss; \
      oss << "(" << comm.getRank() << ") " << ___FILE___ << ":" << __LINE__ << " " << (num) << std::endl;\
      throw(std::bad_alloc(oss.str()); \
    } \
  } \
}

#endif

} //namespace Z2

#endif


