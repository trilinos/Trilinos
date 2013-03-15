
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
#include <Teuchos_RCP.hpp>
#include <Kokkos_DefaultKernels.hpp>
#include <Kokkos_DefaultNode.hpp>
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_ConfigDefs.hpp"

//! Tpetra_CrsMatrixTransposer: A class for transposing an Tpetra_CrsMatrix object.

namespace Tpetra {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
class Map;

/*! \class RowMatrixTransposer
    \brief Construct and (optionally) redistribute the transpose of a CrsMatrix.

    This class is based on the EpetraExt version.  It first transposes the matrix to an
    intermediate version with overlapping row map.  That matrix is then converted to
    a final version whose row map is "unique", i.e., a row is wholly owned by one process.

    This class takes the same template parameters (with the same
    default values) as CrsMatrix.
*/
template <class Scalar, 
	  class LocalOrdinal=int, 
	  class GlobalOrdinal=LocalOrdinal, 
	  class Node=Kokkos::DefaultNode::DefaultNodeType, 
	  class SpMatOps=typename Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps>
class RowMatrixTransposer {
public:
  //! @name Typedefs
  //@{ 
  typedef Scalar scalar_type;
  typedef LocalOrdinal local_ordinal_type;
  typedef GlobalOrdinal global_ordinal_type;
  typedef Node node_type;
  // These match the two typedefs in CrsMatrix.
  typedef SpMatOps mat_vec_type;
  typedef SpMatOps mat_solve_type;

  typedef Map<LocalOrdinal, GlobalOrdinal, Node> map_type;
  typedef CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps> crs_matrix_type;

  //@}
  //! @name Constructor and destructor
  //@{ 

  //! Constructor that takes the matrix to transpose.
  RowMatrixTransposer (const Teuchos::RCP<const crs_matrix_type>& origMatrix);

  /// Constructor that takes the matrix to transpose.
  ///
  /// This method is DEPRECATED, because it is not memory safe.  (If
  /// origMatrix falls out of scope, its reference will be
  /// invalidated.)  Please call the version of the constructor that
  /// takes an <tt>RCP<const crs_matrix_type></tt>.
  TEUCHOS_DEPRECATED RowMatrixTransposer (const crs_matrix_type& origMatrix);

  //@}
  //! @name Forward transformation methods
  //@{ 

  /// Compute and return the transpose of the matrix given to the constructor.
  ///
  /// \param optimizeTranspose [in] If true, optimize the storage of
  ///   the transpose matrix to return.
  ///
  /// \param transposeMatrix [in/out] The target of the transpose
  ///   operation; the matrix into which the result of the transpose
  ///   will be put.
  ///
  /// \param TransposeRowMap [in] If this argument is not null, then
  ///   the transpose matrix will be distributed using this Map as the
  ///   row Map.  If null, this method will evenly distribute the rows
  ///   of the transpose matrix.
  Teuchos::RCP<crs_matrix_type>
  createTranspose (const OptimizeOption optimizeTranspose=DoOptimizeStorage,
		   Teuchos::RCP<const map_type> transposeRowMap=Teuchos::null);
        
private: 
  //! The original matrix to be transposed.
  Teuchos::RCP<const crs_matrix_type> origMatrix_;
};


}

#endif /* TPETRA_ROWMATRIXTRANSPOSER_DECL_HPP */
