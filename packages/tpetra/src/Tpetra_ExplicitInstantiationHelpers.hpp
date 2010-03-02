// @HEADER
// ***********************************************************************
// 
//                    Tpetra: Templated Petra Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
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

#ifdef HAVE_KOKKOS_THRUST
#include <Kokkos_ThrustGPUNode.hpp>
#  if defined(HAVE_KOKKOS_CUDA_FLOAT) && defined(HAVE_TPETRA_INST_FLOAT)
#    define TPETRA_MACRO_TEMPLATE_INSTANT_THRUST_FLOAT(LO,GO,INSTANT_MACRO)\
       INSTANT_MACRO(float,LO,GO,Kokkos::ThrustGPUNode)
#  else
#    define TPETRA_MACRO_TEMPLATE_INSTANT_THRUST_FLOAT(LO,GO,INSTANT_MACRO)
#  endif
#  if defined(HAVE_KOKKOS_CUDA_DOUBLE) && defined(HAVE_TPETRA_INST_DOUBLE)
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
#ifdef HAVE_KOKKOS_TBB
#  include <Kokkos_TBBNode.hpp>
#  define TPETRA_MACRO_TEMPLATE_INSTANT_TBBNODE(SCALAR,LO,GO,INSTANT_MACRO)\
     INSTANT_MACRO(SCALAR,LO,GO,Kokkos::TBBNode)
#else
#  define TPETRA_MACRO_TEMPLATE_INSTANT_TBBNODE(SCALAR,LO,GO,INSTANT_MACRO)
#endif

/** \brief Instantiate a macro template for the Kokkos::TPINode */
#ifdef HAVE_KOKKOS_THREADPOOL
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
