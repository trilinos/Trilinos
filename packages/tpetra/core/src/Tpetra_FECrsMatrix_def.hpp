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
            const Teuchos::RCP<Teuchos::ParameterList>& params) : crs_matrix_type(graph, params)
{
    #ifdef HAVE_TPETRA_DEBUG
    // TODO -->
    // In Tpetra Debug mode, we should check that the graph and offRank graph have no
    // common rows on any given rank (i.e., no intersections).
    #endif

    // Create an offRankMatrix_ object if the offRankGraph isn't null (i.e., if we're not SERIAL)
    if(!offRankGraph.is_null())
    {
        offRankMatrix_ = Teuchos::rcp(new crs_matrix_type(offRankGraph, params));
    }
}


#if 0  // disabled for now
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
FECrsMatrix(const FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& rhs)
{

}
#endif


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
~FECrsMatrix()
{
    // Probably nothing.
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
LocalOrdinal
FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
sumIntoGlobalValues(const GlobalOrdinal globalRow,
                    const Teuchos::ArrayView<const GlobalOrdinal>& cols,
                    const Teuchos::ArrayView<const Scalar>& vals,
                    const bool atomic)
{
    const LocalOrdinal LO_INVALID = Teuchos::OrdinalTraits<LocalOrdinal>::invalid();

    // If globalRow is in 'this' matrix (i.e., is local)
    if(LO_INVALID != this->getRowMap()->getLocalElement(globalRow) )
    {
        this->sumIntoGlobalValues(globalRow, cols, vals, atomic);
    }
    else
    {
        // If not local, let's check the offRankMatrix
        if(LO_INVALID != this->offRankMatrix_->getRowMap()->getLocalElement(globalRow))
        {
            this->offRankMatrix_->sumIntoGlobalValues(globalRow, cols, vals, atomic);
        }
        // If neither local nor offRank, throw an error.
        else
        {
            throw std::runtime_error("FECrsMatrix: sumIntoGlobalValues, globalRow not found.");
        }
    }

  return 0;
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
  return 0;
}


#if 0   // disabled for now
template<class GlobalIndicesViewType, class ImplScalarViewType>
LocalOrdinal
FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
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
    // Call setAllToScalar on the local matrix (this).
    this->setAllToScalar(alpha);

    // Call setAllToSalar on the offRankMatrix if it's not null (i.e., we're not serial)
    if(!offRankMatrix_.is_null())
    {
        offRankMatrix_->setAllToScalar(alpha);
    }
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
globalAssemble()
{
    // Here be dragons...
    // - Will need Tim's new TAFC for this, or equivalent.
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
resumeFill(const Teuchos::RCP<Teuchos::ParameterList>& params)
{
    // Call resumeFill on the local matrix (this).
    this->resumeFill(params);

    // Call resumeFill on the offRankMatrix if it's not null (i.e., we're not serial)
    if(!offRankMatrix_.is_null())
    {
        offRankMatrix_->resumeFill(params);
    }
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>&
FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
operator=(const FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& rhs)
{
  return *this;
}



}  // end namespace Tpetra


//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//
#define TPETRA_FECRSMATRIX_INSTANT(SCALAR,LO,GO,NODE) \
  template class FECrsMatrix<SCALAR, LO, GO, NODE>;



#endif // TPETRA_FECRSMATRIX_DEF
