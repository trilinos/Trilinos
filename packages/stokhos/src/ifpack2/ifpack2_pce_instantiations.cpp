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

#include "Stokhos_Ifpack2.hpp"

typedef Stokhos::StandardStorage<int,double> Storage;
typedef Sacado::PCE::OrthogPoly<double,Storage> pce_type;
typedef pce_type Scalar;
typedef int LocalOrdinal;
typedef int GlobalOrdinal;

#ifdef HAVE_IFPACK2_EXPLICIT_INSTANTIATION
#include "Ifpack2_RILUK_def.hpp"
#include "Ifpack2_ILUT_def.hpp"
#include "Ifpack2_Chebyshev_def.hpp"
#include "Ifpack2_Diagonal_def.hpp"
#include "Ifpack2_Relaxation_def.hpp"
#include "Ifpack2_Krylov_def.hpp"
#ifdef HAVE_TPETRA_EXPLICIT_INSTANTIATION
#include "Tpetra_CrsGraph_def.hpp"
#include "Tpetra_CrsMatrix_def.hpp"
#endif

#include "Ifpack2_ExplicitInstantiationHelpers.hpp"

// Needed to copy these from Ifpack2_Krylov.cpp since they are defined in the
// cpp file
#define IFPACK2_INST_SPARSE_RELAX(S,LO,GO) \
  template class Krylov<Tpetra::CrsMatrix<S,LO,GO>, \
                        Ifpack2::Relaxation<Tpetra::CrsMatrix<S,LO,LO> > >;

#define IFPACK2_INST_SPARSE_ILUT(S,LO,GO) \
  template class Krylov<Tpetra::CrsMatrix<S,LO,GO,>, \
                        Ifpack2::ILUT<Tpetra::CrsMatrix<S,LO,LO> > >;

#define IFPACK2_INST_SPARSE_CHEBY(S,LO,GO) \
  template class Krylov<Tpetra::CrsMatrix<S,LO,GO,>, \
                        Ifpack2::Chebyshev<Tpetra::CrsMatrix<S,LO,LO> > >;

namespace Ifpack2 {
IFPACK2_INST(RILUK,Scalar,LocalOrdinal,GlobalOrdinal)
IFPACK2_INST(ILUT,Scalar,LocalOrdinal,GlobalOrdinal)
IFPACK2_INST(Chebyshev,Scalar,LocalOrdinal,GlobalOrdinal)
IFPACK2_INST(Diagonal,Scalar,LocalOrdinal,GlobalOrdinal)
IFPACK2_INST(Relaxation,Scalar,LocalOrdinal,GlobalOrdinal)
IFPACK2_INST_SPARSE_RELAX(Scalar,LocalOrdinal,GlobalOrdinal)
IFPACK2_INST_SPARSE_ILUT(Scalar,LocalOrdinal,GlobalOrdinal)
IFPACK2_INST_SPARSE_CHEBY(Scalar,LocalOrdinal,GlobalOrdinal)
}

#endif
