// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER

#ifndef TPETRA_EXPLICIT_INSTANTIATION_HELPERS_HPP
#define TPETRA_EXPLICIT_INSTANTIATION_HELPERS_HPP

/*! \file Tpetra_ExplicitInstantiationHelpers.hpp

\brief Macros for helping to explicitly instantiate templated objects.
*/


#include "Tpetra_ConfigDefs.hpp"

//
// 2007/07/10: rabartl: NOTE: Semicolons must only be used at the lowest level
// of final code to ensure that there will not be any empty semicolon lines
// that might issue a compiler warning or error. In general, I like to define
// macros that need a semicolon when you use them because my emacs mode will
// then do the indentation correctly.  However, this is not a big deal since
// these macros only get used in a final *.cpp file and at that point they are
// only used once in the entire mostly empty file.
//

#ifdef HAVE_KOKKOSCLASSIC_THRUST
#include <Kokkos_ThrustGPUNode.hpp>
#  if defined(HAVE_KOKKOSCLASSIC_CUDA_FLOAT) && defined(HAVE_TPETRA_INST_FLOAT)
#    define TPETRA_MACRO_TEMPLATE_INSTANT_THRUST_FLOAT(LO,GO,INSTANT_MACRO)\
       INSTANT_MACRO(float,LO,GO,Kokkos::ThrustGPUNode)
#  else
#    define TPETRA_MACRO_TEMPLATE_INSTANT_THRUST_FLOAT(LO,GO,INSTANT_MACRO)
#  endif
#  if defined(HAVE_KOKKOSCLASSIC_CUDA_DOUBLE) && defined(HAVE_TPETRA_INST_DOUBLE)
#    define TPETRA_MACRO_TEMPLATE_INSTANT_THRUST_DOUBLE(LO,GO,INSTANT_MACRO)\
       INSTANT_MACRO(double,LO,GO,Kokkos::ThrustGPUNode)
#  else
#    define TPETRA_MACRO_TEMPLATE_INSTANT_THRUST_DOUBLE(LO,GO,INSTANT_MACRO)
#  endif
#  define TPETRA_MACRO_TEMPLATE_INSTANT_THRUST_INT(LO,GO,INSTANT_MACRO)\
     INSTANT_MACRO(int,LO,GO,Kokkos::ThrustGPUNode)
#else
#  define TPETRA_MACRO_TEMPLATE_INSTANT_THRUST_FLOAT(LO,GO,INSTANT_MACRO)
#  define TPETRA_MACRO_TEMPLATE_INSTANT_THRUST_DOUBLE(LO,GO,INSTANT_MACRO)
#  define TPETRA_MACRO_TEMPLATE_INSTANT_THRUST_INT(LO,GO,INSTANT_MACRO)
#endif

#define TPETRA_MACRO_TEMPLATE_INSTANT_ALL_SCALAR_THRUSTNODE(LO,GO,MACRO) \
  TPETRA_MACRO_TEMPLATE_INSTANT_THRUST_FLOAT(LO,GO,MACRO) \
  TPETRA_MACRO_TEMPLATE_INSTANT_THRUST_DOUBLE(LO,GO,MACRO) \
  TPETRA_MACRO_TEMPLATE_INSTANT_THRUST_INT(LO,GO,MACRO)

/** \brief Instantiate a macro template for the Kokkos::TBBNode */
#ifdef HAVE_KOKKOSCLASSIC_TBB
#  include <Kokkos_TBBNode.hpp>
#  define TPETRA_MACRO_TEMPLATE_INSTANT_TBBNODE(SCALAR,LO,GO,INSTANT_MACRO)\
     INSTANT_MACRO(SCALAR,LO,GO,Kokkos::TBBNode)
#else
#  define TPETRA_MACRO_TEMPLATE_INSTANT_TBBNODE(SCALAR,LO,GO,INSTANT_MACRO)
#endif

/** \brief Instantiate a macro template for the Kokkos::TPINode */
#ifdef HAVE_KOKKOSCLASSIC_THREADPOOL
#  include <Kokkos_TPINode.hpp>
#  define TPETRA_MACRO_TEMPLATE_INSTANT_TPINODE(SCALAR,LO,GO,INSTANT_MACRO)\
     INSTANT_MACRO(SCALAR,LO,GO,Kokkos::TPINode)
#else
#  define TPETRA_MACRO_TEMPLATE_INSTANT_TPINODE(SCALAR,LO,GO,INSTANT_MACRO)
#endif

/** \brief Instantiate a macro template for the Kokkos::SerialNode */
#include <Kokkos_SerialNode.hpp>
#define TPETRA_MACRO_TEMPLATE_INSTANT_SERIALNODE(SCALAR,LO,GO,INSTANT_MACRO)\
   INSTANT_MACRO(SCALAR,LO,GO,Kokkos::SerialNode)

#define TPETRA_MACRO_TEMPLATE_INSTANT_ALL_CPUNODE(SCALAR,LO,GO,MACRO)\
  TPETRA_MACRO_TEMPLATE_INSTANT_SERIALNODE(SCALAR,LO,GO,MACRO)\
  TPETRA_MACRO_TEMPLATE_INSTANT_TPINODE(SCALAR,LO,GO,MACRO)\
  TPETRA_MACRO_TEMPLATE_INSTANT_TBBNODE(SCALAR,LO,GO,MACRO)

#ifdef HAVE_TPETRA_INST_FLOAT
# define TPETRA_MACRO_TEMPLATE_INSTANT_FLOAT_ALL_CPUNODE(LO,GO,MACRO)\
    TPETRA_MACRO_TEMPLATE_INSTANT_ALL_CPUNODE(float,LO,GO,MACRO)
#else
# define TPETRA_MACRO_TEMPLATE_INSTANT_FLOAT_ALL_CPUNODE(LO,GO,MACRO)
#endif

#ifdef HAVE_TPETRA_INST_DOUBLE
# define TPETRA_MACRO_TEMPLATE_INSTANT_DOUBLE_ALL_CPUNODE(LO,GO,MACRO)\
    TPETRA_MACRO_TEMPLATE_INSTANT_ALL_CPUNODE(double,LO,GO,MACRO)
#else
# define TPETRA_MACRO_TEMPLATE_INSTANT_DOUBLE_ALL_CPUNODE(LO,GO,MACRO)
#endif

#ifdef HAVE_INST_TPETRA_COMPLEX_FLOAT
# define TPETRA_MACRO_TEMPLATE_INSTANT_COMPLEX_FLOAT_ALL_CPUNODE(LO,GO,MACRO)\
    TPETRA_MACRO_TEMPLATE_INSTANT_ALL_CPUNODE(std::complex<float>,LO,GO,MACRO)
#else
# define TPETRA_MACRO_TEMPLATE_INSTANT_COMPLEX_FLOAT_ALL_CPUNODE(LO,GO,MACRO)
#endif

#ifdef HAVE_INST_TPETRA_COMPLEX_DOUBLE
# define TPETRA_MACRO_TEMPLATE_INSTANT_COMPLEX_DOUBLE_ALL_CPUNODE(LO,GO,MACRO)\
    TPETRA_MACRO_TEMPLATE_INSTANT_ALL_CPUNODE(std::complex<double>,LO,GO,MACRO)
#else
# define TPETRA_MACRO_TEMPLATE_INSTANT_COMPLEX_DOUBLE_ALL_CPUNODE(LO,GO,MACRO)
#endif

#define TPETRA_MACRO_TEMPLATE_INSTANT_ALL_SCALAR_CPUNODE(LO,GO,MACRO)\
  TPETRA_MACRO_TEMPLATE_INSTANT_FLOAT_ALL_CPUNODE(LO,GO,MACRO)\
  TPETRA_MACRO_TEMPLATE_INSTANT_DOUBLE_ALL_CPUNODE(LO,GO,MACRO)\
  TPETRA_MACRO_TEMPLATE_INSTANT_COMPLEX_FLOAT_ALL_CPUNODE(LO,GO,MACRO)\
  TPETRA_MACRO_TEMPLATE_INSTANT_COMPLEX_DOUBLE_ALL_CPUNODE(LO,GO,MACRO)\
  TPETRA_MACRO_TEMPLATE_INSTANT_ALL_CPUNODE(int,LO,GO,MACRO)

/** \brief Instantiate a macro template for all Nodes and supported scalar types. */
#define TPETRA_MACRO_TEMPLATE_INSTANT_ALL_SCALAR_NODE(LO,GO,MACRO)\
  TPETRA_MACRO_TEMPLATE_INSTANT_ALL_SCALAR_THRUSTNODE(LO,GO,MACRO)\
  TPETRA_MACRO_TEMPLATE_INSTANT_ALL_SCALAR_CPUNODE(LO,GO,MACRO)

/** \brief Instantiate a macro template for the set of supported scalar, ordinal and node types.
 */
#define TPETRA_MACRO_TEMPLATE_INSTANT_ALL_SCALAR_ORDINAL_NODES(MACRO)\
  TPETRA_MACRO_TEMPLATE_INSTANT_ALL_SCALAR_NODE(short,int,MACRO)\
  TPETRA_MACRO_TEMPLATE_INSTANT_ALL_SCALAR_NODE(int,int,MACRO)

#endif  // TPETRA_EXPLICIT_INSTANTIATION_HELPERS_HPP
