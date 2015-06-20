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

#ifndef XPETRA_CLONER_HPP
#define XPETRA_CLONER_HPP

#include <Kokkos_DefaultNode.hpp>

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_Exceptions.hpp"

#include "Xpetra_Map.hpp"
#ifdef HAVE_XPETRA_TPETRA
#include "Xpetra_TpetraMap.hpp"
#endif

#include "Xpetra_Matrix.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"
#include "Xpetra_MultiVector.hpp"
/** \file Xpetra_Cloner.hpp

Declarations for the class Xpetra::Matrix.
*/
namespace Xpetra {

  template <class LocalOrdinal, class GlobalOrdinal, class Node1, class Node2>
  RCP<Map<LocalOrdinal,GlobalOrdinal,Node2> > clone(const Map<LocalOrdinal,GlobalOrdinal,Node1>& map, const RCP<Node2>& node2) {
    if (map.lib() == UseEpetra)
      throw std::invalid_argument("Map::clone() functionality is only available for Tpetra");
#ifdef HAVE_XPETRA_TPETRA
    RCP<const TpetraMap<LocalOrdinal,GlobalOrdinal,Node1> > tMap = Teuchos::rcp_dynamic_cast<const TpetraMap<LocalOrdinal,GlobalOrdinal,Node1> >(rcpFromRef(map));

    return tMap->clone(node2);
#else
    return Teuchos::null;
#endif
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node1, class Node2>
  RCP<Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node2> > clone(const Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node1>& matrix, const RCP<Node2>& node2) {
    if (matrix.getRowMap()->lib() == UseEpetra)
      throw std::invalid_argument("Matrix::clone() functionality is only available for Tpetra");

    RCP<const CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node1> > tMatrix =
        Teuchos::rcp_dynamic_cast<const CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node1> >(rcpFromRef(matrix));

    return tMatrix->clone(node2);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node1, class Node2>
  RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node2> > clone(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node1>& MV, const RCP<Node2>& node2) {
    if (MV.getMap()->lib() == UseEpetra)
      throw std::invalid_argument("MultiVector::clone() functionality is only available for Tpetra");
#ifdef HAVE_XPETRA_TPETRA
    RCP<const TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node1> > tMV =
        Teuchos::rcp_dynamic_cast<const TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node1> >(rcpFromRef(MV));

    return tMV->clone(node2);
#else
    return Teuchos::null;
#endif
  }


} //namespace Xpetra

#define XPETRA_CLONER_SHORT
#endif //XPETRA_CLONER_HPP
