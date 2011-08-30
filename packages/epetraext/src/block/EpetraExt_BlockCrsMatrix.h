//@HEADER
/*
************************************************************************

              EpetraExt: Extended Linear Algebra Services Package 
                Copyright (2001) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.
 
This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Michael A. Heroux (maherou@sandia.gov) 

************************************************************************
*/
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
