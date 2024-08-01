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

#ifndef TPETRA_MATRIXMATRIX_FWD_HPP
#define TPETRA_MATRIXMATRIX_FWD_HPP

/// \file TpetraExt_MatrixMatrix_fwd.hpp
/// \brief Forward declaration of some Tpetra Matrix Matrix objects
///
/// These are needed by the CrsMatrix to make them friends so they can access CrsMatrix internals

#include "Tpetra_BlockCrsMatrix_fwd.hpp"
#include "Tpetra_CrsMatrix_fwd.hpp"
#include "Tpetra_Import_fwd.hpp"
#include "Tpetra_Map_fwd.hpp"

#include "Teuchos_RCP.hpp"


#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Tpetra {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node> class CrsMatrixStruct;
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node> class BlockCrsMatrixStruct;

namespace MMdetails{

  // Other functions
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
void import_and_extract_views(
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& M,
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > targetMap,
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Mview,
  Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal, Node> > prototypeImporter = Teuchos::null,
  bool userAssertsThereAreNoRemotes = false,
  const std::string& label = std::string(),
  const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
void import_and_extract_views(
  const BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& M,
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > targetMap,
  BlockCrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Mview,
  Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal, Node> > prototypeImporter = Teuchos::null,
  bool userAssertsThereAreNoRemotes = false);

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalOrdinalViewType> struct KernelWrappers;
  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalOrdinalViewType> struct KernelWrappers2;

} // namespace MMdetails
} // namespace Tpetra
#endif // DOXYGEN_SHOULD_SKIP_THIS

#endif // TPETRA_MATRIXMATRIX_FWD_HPP

