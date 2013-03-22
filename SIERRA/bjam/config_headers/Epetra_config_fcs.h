/*
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
//@HEADER
*/

/* $Header$ */

/* Copy of the following file, edited to be used with CMake */
/* src/Epetra_config.h.in.  Generated from configure.ac by autoheader.  */

/* Define the Fortran name mangling to be used for the BLAS */

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef F77_DUMMY_MAIN */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */


/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Define to 1 if you have the <assert.h> header file. */
/* #undef HAVE_ASSERT_H */

/* Define if you have a BLAS library. */
/* #undef HAVE_BLAS */

#define EPETRA_ADDRESS64BIT

/* Define if want to build epetra-abc */
/* #undef HAVE_EPETRA_ARRAY_BOUNDS_CHECK */

/* Define if want to build epetra-examples */
/* #undef HAVE_EPETRA_EXAMPLES */

/* Define if want to build epetra-tests */
/* #undef HAVE_EPETRA_TESTS */

/* Define if want to build with Teuchos enabled */
/* #undef HAVE_EPETRA_TEUCHOS */

/* Define if want to build examples */
/* #undef HAVE_EXAMPLES */

/* Define if you want to build export makefiles. */
/* #undef HAVE_EXPORT_MAKEFILES */

/* Define if want to build with fatal_messages enabled */
/* #undef HAVE_FATAL_MESSAGES */

/* Define if want to build with format_io enabled */
/* #undef HAVE_FORMAT_IO */

/* Define if want to build with Fortran enabled */
#define HAVE_FORTRAN_SUPPORT

/* Define if you are using gnumake - this will shorten your link lines. */
/* #undef HAVE_GNUMAKE */

/* Define if you have LAPACK library. */
/* #undef HAVE_LAPACK */

/* Define if want to build libcheck */
/* #undef HAVE_LIBCHECK */

/* Define to 1 if you have the <math.h> header file. */
/* #undef HAVE_MATH_H */

/* Define to 1 if you have the <memory.h> header file. */
/* #undef HAVE_MEMORY_H */

/* define if we want to use MPI */
#define HAVE_MPI

/* define if we want to use OpenMP */
/* #undef EPETRA_HAVE_OMP */

/* Define if want to build with OSKI enabled */
/* #undef HAVE_OSKI */

/* Define if want to build tests */
/* #undef HAVE_TESTS */

/* Define if want to build with threads enabled */
/* #undef HAVE_THREADS */

/* Define if want to build with warning_messages enabled */
/* #undef HAVE_WARNING_MESSAGES */

/* Define to 1 if you have the ANSI C header files. */
/* #undef STDC_HEADERS */

/* #undef Epetra_ENABLE_CASK */

/* Define if you want to have long long (64 bit) global indices only. */
/* #undef EPETRA_NO_32BIT_GLOBAL_INDICES */

/* Define if you want to have int (32 bit) global indices only. */
/* #undef EPETRA_NO_64BIT_GLOBAL_INDICES */

#ifndef EPETRA_DEPRECATED
#  if (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
#    define EPETRA_DEPRECATED  __attribute__((__deprecated__))
#  else
#    define EPETRA_DEPRECATED
#  endif
#endif

