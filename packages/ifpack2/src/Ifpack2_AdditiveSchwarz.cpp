/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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

#include "Ifpack2_AdditiveSchwarz_decl.hpp"
#include "Ifpack2_ILUT_decl.hpp"
#include "Ifpack2_Details_DenseSolver_decl.hpp"

#if defined(HAVE_IFPACK2_EXPERIMENTAL) && defined(HAVE_IFPACK2_SUPPORTGRAPH)
#include "Ifpack2_SupportGraph_decl.hpp"
#endif

#ifdef HAVE_IFPACK2_EXPLICIT_INSTANTIATION

#include "Ifpack2_AdditiveSchwarz_def.hpp"
#include "Ifpack2_ILUT_def.hpp"
#include "Ifpack2_Details_DenseSolver_def.hpp"

// #if defined(HAVE_IFPACK2_EXPERIMENTAL) && defined(HAVE_IFPACK2_SUPPORTGRAPH)
// #  include "Ifpack2_SupportGraph_def.hpp"
// #endif
#include "Ifpack2_ETIHelperMacros.h"

// mfh 06 Jan 2014: AdditiveSchwarz's second template parameter, the
// type of the subdomain solver, is being deprecated.  It's possible
// already now to control the subdomain solver's type entirely at run
// time, either by specifying it in the input ParameterList, or by
// calling setInnerPreconditioner() with the subdomain solver to use
// (that implements Preconditioner).  Thus, we only need to do
// explicit instantiation for the version of AdditiveSchwarz with one
// template parameter.  For backwards compatibility, we retain
// explicit instantiation for the version with two template
// parameters, only for the commonly used case of ILUT as the
// subdomain solver's type.

// FIXME (mfh 16 Sep 2014) We should really only use RowMatrix here!
// There's no need to instantiate for CrsMatrix too.  All Ifpack2
// preconditioners can and should do dynamic casts if they need a type
// more specific than RowMatrix.

#define IFPACK2_INST_ADDITIVE_SCHWARZ(S,LO,GO) \
  template class AdditiveSchwarz<Tpetra::RowMatrix< S, LO, GO > >; \
  template class AdditiveSchwarz<Tpetra::CrsMatrix< S, LO, GO > >;

#define IFPACK2_INST_ADDITIVE_SCHWARZ_ILUT(S,LO,GO) \
  template class AdditiveSchwarz<Tpetra::RowMatrix< S, LO, GO >, \
                                 Ifpack2::ILUT<Tpetra::RowMatrix< S, LO, GO > > >; \
  template class AdditiveSchwarz<Tpetra::CrsMatrix< S, LO, GO >, \
                                 Ifpack2::ILUT<Tpetra::CrsMatrix< S, LO, GO > > >;


namespace Ifpack2 {

  IFPACK2_ETI_MANGLING_TYPEDEFS()

  IFPACK2_INSTANTIATE_SLG( IFPACK2_INST_ADDITIVE_SCHWARZ )

  IFPACK2_INSTANTIATE_SLG( IFPACK2_INST_ADDITIVE_SCHWARZ_ILUT )

#if defined(HAVE_IFPACK2_KOKKOSCLASSIC) && defined(HAVE_TPETRA_INST_DOUBLE)
#if defined(HAVE_KOKKOSCLASSIC_THRUST) && ! defined(HAVE_KOKKOSCLASSIC_DEFAULTNODE_THRUSTGPUNODE) && defined(HAVE_KOKKOSCLASSIC_CUDA_DOUBLE)
  template class AdditiveSchwarz<Tpetra::CrsMatrix<double, int, int, KokkosClassic::ThrustGPUNode> >;
  template class AdditiveSchwarz<Tpetra::CrsMatrix<double, int, int, KokkosClassic::ThrustGPUNode>,
                                 Ifpack2::ILUT<Tpetra::CrsMatrix<double, int, int, KokkosClassic::ThrustGPUNode> > >;
  template class AdditiveSchwarz<Tpetra::RowMatrix<double, int, int, KokkosClassic::ThrustGPUNode> >;
  template class AdditiveSchwarz<Tpetra::RowMatrix<double, int, int, KokkosClassic::ThrustGPUNode>,
                                 Ifpack2::ILUT<Tpetra::RowMatrix<double, int, int, KokkosClassic::ThrustGPUNode> > >;
#endif

#if defined(HAVE_KOKKOSCLASSIC_THREADPOOL) && ! defined(HAVE_KOKKOSCLASSIC_DEFAULTNODE_TPINODE)
  template class AdditiveSchwarz<Tpetra::CrsMatrix<double, int, int, KokkosClassic::TPINode> >;
  template class AdditiveSchwarz<Tpetra::CrsMatrix<double, int, int, KokkosClassic::TPINode>,
                                 Ifpack2::ILUT<Tpetra::CrsMatrix<double, int, int, KokkosClassic::TPINode> > >;
  template class AdditiveSchwarz<Tpetra::RowMatrix<double, int, int, KokkosClassic::TPINode> >;
  template class AdditiveSchwarz<Tpetra::RowMatrix<double, int, int, KokkosClassic::TPINode>,
                                 Ifpack2::ILUT<Tpetra::RowMatrix<double, int, int, KokkosClassic::TPINode> > >;
#endif
#endif // HAVE_IFPACK2_KOKKOSCLASSIC && HAVE_TPETRA_INST_DOUBLE

}

#endif

