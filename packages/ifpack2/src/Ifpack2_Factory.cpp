/*
//@HEADER
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

#include "Ifpack2_Factory_decl.hpp"

#ifdef HAVE_IFPACK2_EXPLICIT_INSTANTIATION
#  include "Ifpack2_Factory_def.hpp"
#  include "Ifpack2_ExplicitInstantiationHelpers.hpp"
#  include "Ifpack2_ETIHelperMacros.h"
#endif // HAVE_IFPACK2_EXPLICIT_INSTANTIATION

namespace Ifpack2 {

bool supportsUnsymmetric(const std::string& prec_type)
{
  bool result = false;
  if (prec_type == "RELAXATION" ||
      prec_type == "CHEBYSHEV"  ||
      prec_type == "DIAGONAL"   ||
      prec_type == "RILUK"      ||
      prec_type == "ILUT"       ||
      prec_type == "SCHWARZ"    ||
      prec_type == "KRYLOV")
  {
    result = true;
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::invalid_argument, "Ifpack2::supportsUnsymmetric: "
      "Unrecognized preconditioner type prec_type = \"" << prec_type
      << "\"");
  }
  return result;
}

#ifdef HAVE_IFPACK2_EXPLICIT_INSTANTIATION

  // mfh 09 Nov 2013: We can't use the usual IFPACK2_* class macro
  // here because Factory is not a templated class; its methods are.
#define LCLINST(S,LO,GO) \
  template<> \
  Teuchos::RCP<Preconditioner<S, LO, GO, Tpetra::CrsMatrix<S, LO, GO>::node_type> > \
  Factory::create<Tpetra::CrsMatrix<S, LO, GO> > (const std::string& precType, \
                                                  const Teuchos::RCP<const Tpetra::CrsMatrix<S, LO, GO> >& matrix); \
  \
  template<> \
  Teuchos::RCP<Preconditioner<S, LO, GO, Tpetra::CrsMatrix<S, LO, GO>::node_type> > \
  Factory::create<Tpetra::CrsMatrix<S, LO, GO> > (const std::string& precType, \
                                                  const Teuchos::RCP<const Tpetra::CrsMatrix<S, LO, GO> >& matrix, \
                                                  const int overlap);

#if defined(HAVE_KOKKOSCLASSIC_THRUST) && defined(HAVE_KOKKOSCLASSIC_CUDA_DOUBLE) && defined(HAVE_TPETRA_INST_DOUBLE)
  template<>
  Teuchos::RCP<Preconditioner<double, int, int, KokkosClassic::ThrustGPUNode> >
  Factory::create<Tpetra::CrsMatrix<double, int, int, KokkosClassic::ThrustGPUNode> > (const std::string& prec_type,
                                                                                       const Teuchos::RCP<const Tpetra::CrsMatrix<double, int, int, KokkosClassic::ThrustGPUNode> >& matrix);

  template<>
  Teuchos::RCP<Preconditioner<double, int, int, KokkosClassic::ThrustGPUNode> >
  Factory::create<Tpetra::CrsMatrix<double, int, int, KokkosClassic::ThrustGPUNode> > (const std::string& prec_type,
                                                                                       const Teuchos::RCP<const Tpetra::CrsMatrix<double, int, int, KokkosClassic::ThrustGPUNode> >& matrix,
                                                                                       const int overlap);
#endif

#if defined(HAVE_KOKKOSCLASSIC_THREADPOOL) && defined(HAVE_TPETRA_INST_DOUBLE)
  template<>
  Teuchos::RCP<Preconditioner<double, int, int, KokkosClassic::TPINode> >
  Factory::create<Tpetra::CrsMatrix<double, int, int, KokkosClassic::TPINode> > (const std::string& prec_type,
                                                                                 const Teuchos::RCP<const Tpetra::CrsMatrix<double, int, int, KokkosClassic::TPINode> >& matrix);

  template<>
  Teuchos::RCP<Preconditioner<double, int, int, KokkosClassic::TPINode> >
  Factory::create<Tpetra::CrsMatrix<double, int, int, KokkosClassic::TPINode> > (const std::string& prec_type,
                   const Teuchos::RCP<const Tpetra::CrsMatrix<double, int, int, KokkosClassic::TPINode> >& matrix,
                   const int overlap);
#endif

  IFPACK2_ETI_MANGLING_TYPEDEFS()

  IFPACK2_INSTANTIATE_SLG( LCLINST )

#endif

} // namespace Ifpack2

