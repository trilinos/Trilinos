
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

#ifndef TPETRA_ROWMATRIXTRANSPOSER_DECL_HPP
#define TPETRA_ROWMATRIXTRANSPOSER_DECL_HPP

/// \file Tpetra_RowMatrixTransposer_decl.hpp
///
/// Declaration of Tpetra::RowMatrixTransposer.

#include "Tpetra_RowMatrixTransposer_fwd.hpp"
#include "Tpetra_CrsMatrix_fwd.hpp"
#include "Tpetra_Map_fwd.hpp"
#include "Teuchos_RCP.hpp"
#include <string>

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Teuchos {
  // Forward declaration of ParameterList
  class ParameterList;
} // namespace Teuchos
#endif // DOXYGEN_SHOULD_SKIP_THIS

namespace Tpetra {

/// \class RowMatrixTransposer
/// \brief Construct and (optionally) redistribute the explicitly
///   stored transpose of a CrsMatrix.
///
/// This class is based on the EpetraExt version.  It first transposes
/// the matrix to an intermediate version with overlapping row map.
/// That matrix is then converted to a final version whose row map is
/// "unique", i.e., a row is wholly owned by one process.
///
/// This class takes the same template parameters as CrsMatrix.
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
class RowMatrixTransposer {
public:
  //! @name Typedefs
  //@{
  typedef Scalar scalar_type;
  typedef LocalOrdinal local_ordinal_type;
  typedef GlobalOrdinal global_ordinal_type;
  typedef Node node_type;

  typedef Map<LocalOrdinal, GlobalOrdinal, Node> map_type;
  typedef CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crs_matrix_type;

  //@}
  //! @name Constructors
  //@{

  //! Constructor that takes the matrix to transpose.
  RowMatrixTransposer (const Teuchos::RCP<const crs_matrix_type>& origMatrix,const std::string & label = std::string());

  //@}
  //! @name Methods for computing the explicit transpose.
  //@{

  //! Compute and return the transpose of the matrix given to the constructor.
  Teuchos::RCP<crs_matrix_type> createTranspose(const Teuchos::RCP<Teuchos::ParameterList> &params=Teuchos::null);

  /// \brief Compute and return the transpose of the matrix given to the constructor.
  ///
  /// In this call, we (potentially) leave the matrix with an
  /// overlapping row Map.  This is a perfectly valid matrix, but
  /// won't work correctly with some routines in Ifpack or Muelu.
  ///
  /// \warning This routine leaves overlapping rows.  Unless you're
  /// sure that's OK, call createTranspose() instead.
  Teuchos::RCP<crs_matrix_type> createTransposeLocal (const Teuchos::RCP<Teuchos::ParameterList> &params=Teuchos::null);

private:
  //! The original matrix to be transposed.
  Teuchos::RCP<const crs_matrix_type> origMatrix_;

  //! Label for timers
  std::string label_;
};

} // namespace Tpetra

#endif /* TPETRA_ROWMATRIXTRANSPOSER_DECL_HPP */
