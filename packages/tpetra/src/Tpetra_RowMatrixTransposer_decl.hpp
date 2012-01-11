
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


/*! This class provides capabilities to construct a transpose matrix of an existing Tpetra_CrsMatrix
	  object and (optionally) redistribute it across a parallel distributed memory machine.
*/
template <class Scalar, 
	class LocalOrdinal=int, 
	class GlobalOrdinal=LocalOrdinal, 
	class Node=Kokkos::DefaultNode::DefaultNodeType, 
    class SpMatOps   = typename Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps >
class RowMatrixTransposer {
    
  public:

  //! @name Constructors/destructors
  //@{ 
  //! Primary Tpetra_CrsMatrixTransposer constructor.
  /*!
    \param origMatrix An existing Tpetra_CrsMatrix object.  The Tpetra_CrsMatrix, the LHS and RHS pointers
		       do not need to be defined before this constructor is called.

    \return Pointer to a Tpetra_CrsMatrixTransposer object.

  */ 
  RowMatrixTransposer(const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& origMatrix);

  //! Tpetra_CrsMatrixTransposer destructor.
  
  virtual ~RowMatrixTransposer();
  //@}
  
  //! @name Forward transformation methods
  //@{ 
  
  //! Generate a new Tpetra_CrsMatrix as the transpose of an Tpetra_CrsMatrix passed into the constructor.
  /*! Constructs a new Tpetra_CrsMatrix that is a copy of the Tpetra_CrsMatrix passed in to the constructor.
		
		\param optimizeTranspose Optimizes the storage of the newly created Transpose matrix
		\param transposeMatrix The matrix in which the result of the tranpose operation will be put.
		\param TransposeRowMap If this argument is defined, the transpose matrix will be distributed
		       using this map as the row map for the transpose.	If null, the function will evenly distribute
			   the rows of the tranpose matrix.
  */
  RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps> >
  createTranspose(const OptimizeOption optimizeTranspose=DoOptimizeStorage
    , Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > transposeRowMap = Teuchos::null);

	
 private: 
	//The original matrix to be transposed.
	const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& origMatrix_;
	//The matrix in which the result of the tranpose is placed.
	RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps> > transposeMatrix_;
	//Whether or not to optimize the storage of the transpose matrix.
	OptimizeOption optimizeTranspose_;	
	const RCP<const Teuchos::Comm<int> > comm_;
	GlobalOrdinal indexBase_;
};


}

#endif /* TPETRA_ROWMATRIXTRANSPOSER_DECL_HPP */
