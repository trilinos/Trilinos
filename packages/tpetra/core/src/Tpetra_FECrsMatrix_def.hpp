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

#ifndef TPETRA_FECRSMATRIX_DEF_HPP
#define TPETRA_FECRSMATRIX_DEF_HPP

#include "Tpetra_CrsMatrix.hpp"



namespace Tpetra {



template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
FECrsMatrix(const Teuchos::RCP<const crs_graph_type>& graph,
            const Teuchos::RCP<const crs_graph_type>& offRankGraph,
            const Teuchos::RCP<Teuchos::ParameterList>& params)
{
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
~FECrsMatrix()
{
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
LocalOrdinal
FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
sumIntoGlobalValues(const GlobalOrdinal globalRow,
                    const Teuchos::ArrayView<const GlobalOrdinal>& cols,
                    const Teuchos::ArrayView<const Scalar>& vals,
                    const bool atomic)
{
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
LocalOrdinal
FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
sumIntoGlobalValues(const GlobalOrdinal globalRow,
                    const LocalOrdinal numEnt,
                    const Scalar vals[],
                    const GlobalOrdinal cols[],
                    const bool atomic)
{
}


#if 0
template<class GlobalIndicesViewType, class ImplScalarViewType>
LocalOrdinal
sumIntoGlobalValues (const GlobalOrdinal globalRow,
                     const typename UnmanagedView<GlobalIndicesViewType>::type& inputInds,
                     const typename UnmanagedView<ImplScalarViewType>::type& inputVals,
                     const bool atomic)
{
}
#endif


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    setAllToScalar(const Scalar &alpha)
{
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
globalAssemble()
{
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
resumeFill(const Teuchos::RCP<Teuchos::ParameterList>& params)
{
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
FECrsMatrix(const FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& rhs)
{
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>&
FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
operator=(const FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& rhs)
{
}



}  // end namespace Tpetra



#endif // TPETRA_FECRSMATRIX_DEF
