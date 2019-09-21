// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

// WARNING: This code is experimental. Backwards compatibility should not be expected.

#ifndef XPETRA_MATRIXVIEW_HPP
#define XPETRA_MATRIXVIEW_HPP

#include <Teuchos_Describable.hpp>
#include <Kokkos_DefaultNode.hpp>

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_Map.hpp"

/** \file Xpetra_MatrixView.hpp

Declarations for the class Xpetra::MatrixView.
*/
namespace Xpetra {

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node = KokkosClassic::DefaultNode::DefaultNodeType>
  class MatrixView { // TODO : virtual public Teuchos::Describable {
    typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;

  public:

    //! @name Constructor/Destructor Methods
    //@{

    //! Constructor
    MatrixView(const RCP<const Map> &rowMap, const RCP<const Map> &colMap)
      : rowMap_ (rowMap), colMap_(colMap), maxEigValueEstimate_(-Teuchos::ScalarTraits<Scalar>::one())
    { }

    //! Destructor
    virtual ~MatrixView() {}

    //@}

    //! @name Map access methods
    //@{
    //! Returns the Map that describes the row distribution in this matrix.
    const RCP<const Map> & GetRowMap() const { return rowMap_; }

    //! \brief Returns the Map that describes the column distribution in this matrix.
    const RCP<const Map> & GetColMap() const { return colMap_; }

    //! Returns the Map that describes the row distribution in this matrix.
    void SetRowMap(const RCP<const Map> & rowMap) { rowMap_ = rowMap; }

    //! \brief Set the Map that describes the column distribution in this matrix.
    void SetColMap(const RCP<const Map> & colMap) { colMap_ = colMap; }
    //@}

    //! \brief Set an maximum eigenvalue estimate for this matrix.
    void SetMaxEigenvalueEstimate(Scalar const &sigma) { maxEigValueEstimate_ = sigma; }

    //! \brief Return the maximum eigenvalue estimate for this matrix.
    Scalar GetMaxEigenvalueEstimate() const { return maxEigValueEstimate_; }

  private:
    RCP<const Map> rowMap_;
    RCP<const Map> colMap_;

    Scalar maxEigValueEstimate_;

  }; // class MatrixView

} // namespace Xpetra

#define XPETRA_MATRIXVIEW_SHORT
#endif //XPETRA_MATRIX_VIEW_DECL_HPP
