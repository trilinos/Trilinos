// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_MUELU_PRECONDITIONER_HPP
#define STOKHOS_MUELU_PRECONDITIONER_HPP

#include "Teuchos_RCP.hpp"
#include "SGPreconditioner.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Tpetra_MultiVector.hpp"
#include "MueLu_CreateTpetraPreconditioner.hpp"

namespace Kokkos {
namespace Example {

  template<class S, class LO, class GO, class N>
  class MueLuPreconditioner :
    public SGPreconditioner<S, LO, GO, N> {

  public:

    //! Constructor
    MueLuPreconditioner() {}

    //! Destructor
    virtual ~MueLuPreconditioner() {}

    //! Setup preconditioner
    virtual
    Teuchos::RCP<Tpetra::Operator<S,LO,GO,N> >
    setupPreconditioner(
      const Teuchos::RCP<Tpetra::CrsMatrix<S,LO,GO,N> >& A,
      const std::string& xmlFileName,
      const Teuchos::RCP<Tpetra::MultiVector<double,LO,GO,N> >& coords);

  private:

    //! Private to prohibit copying
    MueLuPreconditioner(const MueLuPreconditioner&);

    //! Private to prohibit copying
    MueLuPreconditioner& operator=(const MueLuPreconditioner&);


  }; // class MeanBasedPreconditioner

}
}

template<class S, class LO, class GO, class N>
Teuchos::RCP<Tpetra::Operator<S,LO,GO,N> >
Kokkos::Example::MueLuPreconditioner<S,LO,GO,N>::
setupPreconditioner(const Teuchos::RCP<Tpetra::CrsMatrix<S,LO,GO,N> >& A,
                    const std::string& xmlFileName,
                    const Teuchos::RCP<Tpetra::MultiVector<double,LO,GO,N> >& coords)
{
  typedef MueLu::TpetraOperator<S,LO,GO,N> PreconditionerType;
  Teuchos::RCP<PreconditionerType> mueluPreconditioner;
  mueluPreconditioner = MueLu::CreateTpetraPreconditioner(A,xmlFileName,coords);
  return mueluPreconditioner;
}


#endif // STOKHOS_MUELU_PRECONDITIONER_HPP
