/*
//@HEADER
// ************************************************************************
// 
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/


/* Define the Fortran name mangling to be used for the BLAS */
#define F77_BLAS_MANGLE(name,NAME) name ## _

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef F77_DUMMY_MAIN */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define F77_FUNC(name,NAME) name ## _

/* As F77_FUNC, but for C identifiers containing underscores. */
#define F77_FUNC_(name,NAME) name ## _

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Define to 1 if you have the <fpu_control.h> header file. */
/* #undef HAVE_FPU_CONTROL_H */

/* define if the compiler supports abi::__cxa_demangle(...) */
#define HAVE_GCC_ABI_DEMANGLE

/* define if we want to use MPI */
/* #undef HAVE_MPI */

/* #undef HAVE_COMPLEX_BLAS_PROBLEM */

/* #undef HAVE_FIXABLE_COMPLEX_BLAS_PROBLEM */

/* define if the compiler supports access of protected templated nested
   classes in derived classes */
/* #undef HAVE_PROTECTED_NESTED_TEMPLATE_CLASS_ACCESS */

/* Define if want to build teuchos-abc */
/* #undef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK */

/* Define if want to build teuchos-blasfloat */
#define HAVE_TEUCHOS_BLASFLOAT

/* Define if want to build teuchos-boost */
/* #undef HAVE_TEUCHOS_BOOST */

/* Define if you want to build teuchos-qt */
/* #undef HAVE_TEUCHOS_QT */

/* Define if want to build teuchos-qd */
/* #undef HAVE_TEUCHOS_QD */

/* Define if ScalarTraits promotes double to QD */
/* #undef HAVE_TEUCHOS_DOUBLE_TO_QD */

/* Define if want to build teuchos-arprec */
/* #undef HAVE_TEUCHOS_ARPREC */

/* Define if ScalarTraits promotes double to ARPREC */
/* #undef HAVE_TEUCHOS_DOUBLE_TO_ARPREC */

/* Define if want to build teuchos-comm_timers */
/* #undef HAVE_TEUCHOS_COMM_TIMERS */

/* #undef HAVE_TEUCHOS_FLOAT */

#define TEUCHOS_ORDINAL_TYPE ptrdiff_t

/* Define if want to build teuchos-complex */
/* #undef HAVE_TEUCHOS_COMPLEX */

/* Define if want to build teuchos-long-long */
/* #undef HAVE_TEUCHOS_LONG_LONG_INT */

/* Define if want to build teuchos-debug */
/* #undef HAVE_TEUCHOS_DEBUG */

/* #undef HAVE_TEUCHOS_DEBUG_RCP_NODE_TRACING */

/* #undef HAS_TEUCHOS_BOOST_IS_POLYMORPHIC */

/* Define if want to build teuchos-demangle */
#define HAVE_TEUCHOS_DEMANGLE

/* Define if want to build teuchos-expat */
/* #undef HAVE_TEUCHOS_EXPAT */

/* Define if want to build teuchos-explicit_instantiation */
/* #undef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION */

/* Define if want to build teuchos-extended */
#define HAVE_TEUCHOS_EXTENDED

/* Define if want to build teuchos-gmp */
/* #undef HAVE_TEUCHOS_GNU_MP */

/* Define if want to build teuchos-libxml2 */
/* #undef HAVE_TEUCHOS_LIBXML2 */

/* #undef HAVE_TEUCHOS_C_EXCEPTIONS */

/* #undef HAVE_TEUCHOS_LINK */

/* #undef HAVE_TEUCHOS_BFD */

/* #undef HAVE_TEUCHOS_STACKTRACE */

/* #undef HAVE_TEUCHOS_DEFAULT_STACKTRACE */

/* template qualifier required for calling template methods from non-template
   code */
#define INVALID_TEMPLATE_QUALIFIER 

#ifndef TEUCHOS_DEPRECATED
#  if (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
#    define TEUCHOS_DEPRECATED  __attribute__((__deprecated__))
#  else
#    define TEUCHOS_DEPRECATED
#  endif
#endif

