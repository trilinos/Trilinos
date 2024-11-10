/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef MACROS_DH
#define MACROS_DH

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif

#ifndef FMAX
#define FMAX(a,b)  ((FABS(a)) > (FABS(b)) ? (FABS(a)) : (FABS(b)))
#endif

#ifndef MAX
#define MAX(a,b)   ((a) > (b) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a,b) ((a)<(b)?(a):(b))
#endif

#ifndef ABS
#define ABS(x) (((x)<0)?(-(x)):(x))
#endif

#ifndef FABS
#define FABS(a)    ((a) < 0 ? -(a) : a)
#endif

/* used in Mat_SEQ_PrintTriples, so matlab won't discard zeros (yuck!) */
#define _MATLAB_ZERO_  1e-100


/*---------------------------------------------------------------------- 
 * macros for error handling everyplace except in main.
 *---------------------------------------------------------------------- */

/* for future expansion: should check that "ptr" points to
   a valid memory address, if not null.
*/
#define ASSERT_DH(ptr) \
    { \
      if (ptr == NULL) { \
        sprintf(msgBuf_dh, "%s is NULL", ptr); \
        SET_V_ERROR(msgBuf_dh); \
      } \
    }


#if 0
#define CHECK_MPI_V_ERROR(errCode)  \
      { \
        if (errCode) { \
          int len; \
          MPI_Error_string(errCode, msgBuf_dh, &len); \
          setError_dh(msgBuf_dh, __FUNC__, __FILE__, __LINE__); \
          return; \
        } \
      }

#define CHECK_MPI_ERROR(errCode)  \
      { \
        if (errCode) { \
          int len; \
          MPI_Error_string(errCode, msgBuf_dh, &len); \
          setError_dh(msgBuf_dh, __FUNC__, __FILE__, __LINE__); \
          return(errCode); \
        } \
      }
#endif

#define CHECK_MPI_V_ERROR(errCode)  \
      { \
        if (errCode) { \
          setError_dh("MPI error!", __FUNC__, __FILE__, __LINE__); \
          printErrorMsg(stderr); \
          MPI_Abort(comm_dh, -1); \
        } \
      }

#define CHECK_MPI_ERROR(errCode)  \
      { \
        if (errCode) { \
          setError_dh("MPI error!", __FUNC__, __FILE__, __LINE__); \
          printErrorMsg(stderr); \
          MPI_Abort(comm_dh, -1); \
        } \
      }

#define SET_V_ERROR(msg)  \
      { setError_dh(msg, __FUNC__, __FILE__, __LINE__); \
        printErrorMsg(stderr); \
        MPI_Abort(comm_dh, -1); \
      }

#define SET_ERROR(retval, msg)  \
      { setError_dh(msg, __FUNC__, __FILE__, __LINE__); \
        printErrorMsg(stderr); \
        MPI_Abort(comm_dh, -1); \
      }

#define CHECK_V_ERROR   \
          if (errFlag_dh) { \
            setError_dh("",  __FUNC__, __FILE__, __LINE__); \
            printErrorMsg(stderr); \
            MPI_Abort(comm_dh, -1); \
          }

#define CHECK_ERROR(retval)  \
          if (errFlag_dh) { \
            setError_dh("",  __FUNC__, __FILE__, __LINE__); \
            printErrorMsg(stderr); \
            MPI_Abort(comm_dh, -1); \
          }

/*---------------------------------------------------------------------- 
 * informational macros
 *---------------------------------------------------------------------- */

#define SET_INFO(msg)  setInfo_dh(msg, __FUNC__, __FILE__, __LINE__);

/*---------------------------------------------------------------------- 
 * macros for tracking the function call stack
 *---------------------------------------------------------------------- */
#ifdef OPTIMIZED_DH

#define START_FUNC_DH   \
          dh_StartFunc(__FUNC__, __FILE__, __LINE__, 1);  \
          {

#define END_FUNC_DH     \
          } \
          dh_EndFunc(__FUNC__, 1);

#define END_FUNC_VAL(a) \
         dh_EndFunc(__FUNC__, 1); \
         return a ; \
         }

#define START_FUNC_DH_2 /**/
#define END_FUNC_DH_2 /**/
#define END_FUNC_VAL_2(a) return a ;
#else

#define START_FUNC_DH  \
          dh_StartFunc(__FUNC__, __FILE__, __LINE__, 1); \
          if (logFuncsToStderr || logFuncsToFile)\
            Error_dhStartFunc(__FUNC__, __FILE__, __LINE__); \
          {

#define END_FUNC_DH   \
          dh_EndFunc(__FUNC__, 1); \
          if (logFuncsToStderr || logFuncsToFile) \
            Error_dhEndFunc(__FUNC__); \
          return; \
          } \

#define START_FUNC_DH_2  \
          dh_StartFunc(__FUNC__, __FILE__, __LINE__, 2); \
          if (logFuncsToStderr || logFuncsToFile)\
            Error_dhStartFunc(__FUNC__, __FILE__, __LINE__); \
          {

#define END_FUNC_DH_2   \
          dh_EndFunc(__FUNC__, 2); \
          if (logFuncsToStderr || logFuncsToFile) \
            Error_dhEndFunc(__FUNC__); \
          return; \
          } \


#define END_FUNC_VAL(retval) \
          dh_EndFunc(__FUNC__, 1); \
          if (logFuncsToStderr || logFuncsToFile) \
            Error_dhEndFunc(__FUNC__); \
          return(retval); \
          } \

#define END_FUNC_VAL_2(retval) \
          dh_EndFunc(__FUNC__, 2); \
          if (logFuncsToStderr || logFuncsToFile) \
            Error_dhEndFunc(__FUNC__); \
          return(retval); \
          } \


#endif

#endif /* #ifndef MACROS_DH */
