//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
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
// ***********************************************************************
//@HEADER

#ifndef EPETRAEXT_BLOCKCRSMATRIX_H
#define EPETRAEXT_BLOCKCRSMATRIX_H

#include <vector>

#include "Epetra_CrsMatrix.h" 

//! EpetraExt::BlockCrsMatrix: A class for constructing a distributed block matrix.

/*! The EpetraExt::BlockCrsMatrix allows construction of a block matrix made up of Epetra_CrsMatrix blocks as well as access to the full systems as a Epetra_CrsMatrix.  It derives from and extends the Epetra_CrsMatrix class

<b>Constructing EpetraExt::BlockCrsMatrix objects</b>

*/    

namespace EpetraExt {

class BlockCrsMatrix: public Epetra_CrsMatrix {
 public:

  //@{ \name Constructors/Destructor.
  //! BlockCrsMatrix constuctor with one block row per processor.
  /*! Creates a BlockCrsMatrix object and allocates storage.  
    
	\param In
	BaseGraph - Graph determining individual block structure, can be distrib. over subset of proc.'s
	\param In 
	RowStencil - Describes the stencil for block row on this processor (i.e. (-1 0 1) centered difference)
	\param In
	RowIndex - Defines the index used for this block row.
  */
  BlockCrsMatrix( const Epetra_CrsGraph & BaseGraph, const std::vector<int> & RowStencil, int RowIndex, const Epetra_Comm & GlobalComm );
  
  //! BlockCrsMatrix constuctor with multiple block rows per processor.
  /*! Creates a BlockCrsMatrix object and allocates storage.  
    
	\param In
	BaseGraph - Graph determining individual block structure, can be distrib. over subset of proc.'s
	\param In 
	RowStencil - Describes the stencil for block row on this processor (i.e. (-1 0 1) centered difference)
	\param In
	RowIndices - Defines the indices used for this block row.
  */
  BlockCrsMatrix( const Epetra_CrsGraph & BaseGraph, const std::vector< std::vector<int> > & RowStencil, const std::vector<int> & RowIndices, const Epetra_Comm & GlobalComm );

  //! Version taking a local block graph
  BlockCrsMatrix( const Epetra_CrsGraph & BaseGraph, const Epetra_CrsGraph& LocalBlockGraph, const Epetra_Comm & GlobalComm );

  BlockCrsMatrix( const Epetra_RowMatrix & BaseMatrix, const std::vector< std::vector<int> > & RowStencil, const std::vector<int> & RowIndices, const Epetra_Comm & GlobalComm );
  
  //! Copy constructor.
  BlockCrsMatrix( const BlockCrsMatrix & Matrix );

  //! Destructor
  virtual ~BlockCrsMatrix();
  //@}
  
  //! Local Stencil Info
  const std::vector<int> & Stencil( int i = 0 ) { return RowStencil_[i]; }

  //! RowIndex
  int RowIndex( int i = 0 ) { return RowIndices_[i]; }
	
  //! Routine for loading a base matrices values into the large Block Matrix
  //! The Row and Col arguments are indices into RowStencil 
  void LoadBlock(const Epetra_RowMatrix & BaseMatrix, const int Row, const int Col);

  //! Routine for summing base matrices values into the large Block Matrix
  //! The Row and Col arguments are indices into RowStencil 
  void SumIntoBlock(double alpha, const Epetra_RowMatrix & BaseMatrix, const int Row, const int Col);

  //! Routine for summing base matrices values into the large Block Matrix
  //! The Row and Col arguments are global indices
  void SumIntoGlobalBlock(double alpha, const Epetra_RowMatrix & BaseMatrix, const int Row, const int Col);

  //! Sum Entries into Block matrix using base-matrix numbering plus block Row and Col
  //! The Row and Col arguments are indices into RowStencil 
  void BlockSumIntoGlobalValues(const int BaseRow, int NumIndices,
     double* Values, const int* Indices, const int Row, const int Col);
  void BlockReplaceGlobalValues(const int BaseRow, int NumIndices,
     double* Values, const int* Indices, const int Row, const int Col);
  void BlockExtractGlobalRowView(const int BaseRow, int& NumEntries, 
     double*& Values, const int Row, const int Col);

  void ExtractBlock(Epetra_CrsMatrix & BaseMatrix, const int Row, const int Col);

 protected:

  Epetra_CrsGraph BaseGraph_;

  std::vector< std::vector<int> > RowStencil_;

  std::vector<int> RowIndices_; 

  int ROffset_;
  int COffset_;

};

} //namespace EpetraExt

#endif /* EPETRA_CRSMATRIX_H */
