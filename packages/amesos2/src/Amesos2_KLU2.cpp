// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2011 Sandia Corporation
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
// ***********************************************************************
//
// @HEADER

#include "Amesos2_KLU2_decl.hpp"

#ifdef HAVE_AMESOS2_EXPLICIT_INSTANTIATION

#  include "Amesos2_KLU2_def.hpp"
#  include "Amesos2_ExplicitInstantiationHelpers.hpp"

namespace Amesos2 {
#ifdef HAVE_AMESOS2_EPETRA
  AMESOS2_SOLVER_EPETRA_INST(KLU2);
#endif

#ifdef HAVE_TPETRA_INST_FLOAT
  AMESOS2_SOLVER_TPETRA_INST(KLU2,float,int,int);
#endif
#ifdef HAVE_TPETRA_INST_DOUBLE
  AMESOS2_SOLVER_TPETRA_INST(KLU2,double,int,int);
#endif
/*#ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
  AMESOS2_SOLVER_TPETRA_INST(KLU2,std::complex<float>,int,int);
#endif*/
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
  AMESOS2_SOLVER_TPETRA_INST(KLU2,std::complex<double>,int,int);
#endif
#ifdef HAVE_TPETRA_INST_INT_LONG
  AMESOS2_SOLVER_TPETRA_INST(KLU2,double,int,long);
#endif
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
  AMESOS2_SOLVER_TPETRA_INST(KLU2,double,int,long long int);
#endif
}

//
// 7-Nov-2014: JJH code copied from Amesos2_SuperLU.cpp.
//
#include "Kokkos_DefaultNode.hpp"
#include "Tpetra_ETIHelperMacros.h"

#define AMESOS2_KLU2_LOCAL_INSTANT(S,LO,GO,N)                        \
  template class Amesos2::KLU2<Tpetra::CrsMatrix<S, LO, GO, N>,      \
                                  Tpetra::MultiVector<S, LO, GO,  N> >;

TPETRA_ETI_MANGLING_TYPEDEFS()

#if defined(HAVE_KOKKOSCLASSIC_THREADPOOL) && defined(HAVE_TPETRA_INST_DOUBLE) && !defined(HAVE_KOKKOSCLASSIC_DEFAULTNODE_TPINODE)
  AMESOS2_KLU2_LOCAL_INSTANT(double, int, int, KokkosClassic_TPINode)
#endif

#if defined(HAVE_KOKKOSCLASSIC_THRUST) && defined(HAVE_TPETRA_INST_DOUBLE) && !defined(HAVE_KOKKOSCLASSIC_DEFAULTNODE_THRUSTGPUNODE)
  AMESOS2_KLU2_LOCAL_INSTANT(double, int, int, KokkosClassic_ThrustGPUNode)
#endif

#if defined(HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT) && !defined(HAVE_KOKKOSCLASSIC_DEFAULTNODE_SERIALWRAPPERNODE) && defined(HAVE_TPETRA_INST_DOUBLE) && defined(TPETRA_HAVE_KOKKOS_REFACTOR)
  AMESOS2_KLU2_LOCAL_INSTANT(double, int, int, Kokkos_Compat_KokkosSerialWrapperNode)
#endif

#if defined(HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT) && defined(KOKKOS_HAVE_PTHREAD) && !defined(HAVE_KOKKOSCLASSIC_DEFAULTNODE_THREADSWRAPPERNODE) && defined(HAVE_TPETRA_INST_DOUBLE) && defined(TPETRA_HAVE_KOKKOS_REFACTOR)
  AMESOS2_KLU2_LOCAL_INSTANT(double, int, int, Kokkos_Compat_KokkosThreadsWrapperNode)
#endif

#if defined(HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT) && defined(KOKKOS_HAVE_OPENMP) && !defined(HAVE_KOKKOSCLASSIC_DEFAULTNODE_OPENMPWRAPPERNODE) && defined(HAVE_TPETRA_INST_DOUBLE) && defined(TPETRA_HAVE_KOKKOS_REFACTOR)
  AMESOS2_KLU2_LOCAL_INSTANT(double, int, int, Kokkos_Compat_KokkosOpenMPWrapperNode)
#endif

#if defined(HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT) && defined(KOKKOS_HAVE_CUDA) && !defined(HAVE_KOKKOSCLASSIC_DEFAULTNODE_CUDAWRAPPERNODE) && defined(HAVE_TPETRA_INST_DOUBLE) && defined(TPETRA_HAVE_KOKKOS_REFACTOR)
  AMESOS2_KLU2_LOCAL_INSTANT(double, int, int, Kokkos_Compat_KokkosCudaWrapperNode)
#endif

#endif  // HAVE_AMESOS2_EXPLICIT_INSTANTIATION
