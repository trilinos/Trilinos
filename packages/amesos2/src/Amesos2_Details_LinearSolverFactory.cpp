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

// We need both the decl and def files here, whether or not ETI is on.
#include "Amesos2_Details_LinearSolverFactory_decl.hpp"
#include "Amesos2_Details_LinearSolverFactory_def.hpp"
// We need this whether or not ETI is on, in order to define typedefs
// for making Tpetra's macros work.
#include "TpetraCore_ETIHelperMacros.h"

// Define typedefs that make the Tpetra macros work.
TPETRA_ETI_MANGLING_TYPEDEFS()

#ifdef HAVE_AMESOS2_EXPLICIT_INSTANTIATION

#ifdef HAVE_AMESOS2_EPETRA
// Do explicit instantiation of Amesos2::Details::LinearSolverFactory, for Epetra objects.
template class Amesos2::Details::LinearSolverFactory<Epetra_MultiVector, Epetra_Operator>;
#endif // HAVE_AMESOS2_EPETRA

// mfh 23 Jul 2015: Amesos2 has a required dependency on Tpetra,
// so we don't have to protect use of Tpetra with a macro.

// Macro for doing explicit instantiation of
// Amesos2::Details::LinearSolverFactory, for Tpetra objects, with given
// Tpetra template parameters (Scalar, LocalOrdinal, GlobalOrdinal,
// Node).
#define LCLINST(SC, LO, GO, NT) \
  template class Amesos2::Details::LinearSolverFactory<Tpetra::MultiVector<SC, LO, GO, NT>, \
                                                 Tpetra::Operator<SC, LO, GO, NT> >;

// Do explicit instantiation of Amesos2::Details::LinearSolverFactory, for
// Tpetra objects, for all combinations of Tpetra template parameters
// for which Tpetra does explicit template instantiation (ETI).
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( LCLINST )

#endif // HAVE_AMESOS2_EXPLICIT_INSTANTIATION

//
// Register Amesos2's solver factory/ies with the central registry.
// Do this whether or not ETI is on.
//

namespace Amesos2 {
namespace Details {

template<class MV, class OP>
class RegisterLinearSolverFactory {
public:
  RegisterLinearSolverFactory () {
#ifdef HAVE_TEUCHOSCORE_CXX11
    typedef std::shared_ptr<Amesos2::Details::LinearSolverFactory<MV, OP> > ptr_type;
    //typedef std::shared_ptr<Trilinos::Details::LinearSolverFactory<MV, OP> > base_ptr_type;
#else
    typedef Teuchos::RCP<Amesos2::Details::LinearSolverFactory<MV, OP> > ptr_type;
    //typedef Teuchos::RCP<Trilinos::Details::LinearSolverFactory<MV, OP> > base_ptr_type;
#endif // HAVE_TEUCHOSCORE_CXX11

    ptr_type factory (new Amesos2::Details::LinearSolverFactory<MV, OP> ());
    Trilinos::Details::registerLinearSolverFactory<MV, OP> ("Amesos2", factory);
  }
};

} // namespace Details
} // namespace Trilinos

namespace { // (anonymous)

#ifdef HAVE_AMESOS2_EPETRA
  Amesos2::Details::RegisterLinearSolverFactory<Epetra_MultiVector, Epetra_Operator>
  registerer_Epetra;

#endif // HAVE_AMESOS2_EPETRA

#define AMESOS2_DETAILS_REGISTER(SC, LO, GO, NT) \
  Amesos2::Details::RegisterLinearSolverFactory<Tpetra::MultiVector<SC, LO, GO, NT>, \
                                                Tpetra::Operator<SC, LO, GO, NT> > \
    registerer_Tpetra_##SC##_##LO##_##GO##_##NT ;

TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( AMESOS2_DETAILS_REGISTER )

} // namespace (anonymous)
