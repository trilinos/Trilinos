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

#if defined(EpetraExt_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The EpetraExt package is deprecated"
#endif
#endif

#include <vector>

#include "Epetra_ConfigDefs.h" 
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
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  BlockCrsMatrix( const Epetra_CrsGraph & BaseGraph, const std::vector<int> & RowStencil, int RowIndex, const Epetra_Comm & GlobalComm );
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  BlockCrsMatrix( const Epetra_CrsGraph & BaseGraph, const std::vector<long long> & RowStencil, long long RowIndex, const Epetra_Comm & GlobalComm );
#endif

  //! BlockCrsMatrix constuctor with multiple block rows per processor.
  /*! Creates a BlockCrsMatrix object and allocates storage.  
    
	\param In
	BaseGraph - Graph determining individual block structure, can be distrib. over subset of proc.'s
	\param In 
	RowStencil - Describes the stencil for block row on this processor (i.e. (-1 0 1) centered difference)
	\param In
	RowIndices - Defines the indices used for this block row.
  */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  BlockCrsMatrix( const Epetra_CrsGraph & BaseGraph, const std::vector< std::vector<int> > & RowStencil, const std::vector<int> & RowIndices, const Epetra_Comm & GlobalComm );
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  BlockCrsMatrix( const Epetra_CrsGraph & BaseGraph, const std::vector< std::vector<long long> > & RowStencil, const std::vector<long long> & RowIndices, const Epetra_Comm & GlobalComm );
#endif

  //! Version taking a local block graph
  BlockCrsMatrix( const Epetra_CrsGraph & BaseGraph, const Epetra_CrsGraph& LocalBlockGraph, const Epetra_Comm & GlobalComm );

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  BlockCrsMatrix( const Epetra_RowMatrix & BaseMatrix, const std::vector< std::vector<int> > & RowStencil, const std::vector<int> & RowIndices, const Epetra_Comm & GlobalComm );
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  BlockCrsMatrix( const Epetra_RowMatrix & BaseMatrix, const std::vector< std::vector<long long> > & RowStencil, const std::vector<long long> & RowIndices, const Epetra_Comm & GlobalComm );
#endif
  
  //! Copy constructor.
  BlockCrsMatrix( const BlockCrsMatrix & Matrix );

  //! Destructor
  virtual ~BlockCrsMatrix();
  //@}
  
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  //! Local Stencil Info
  const std::vector<int> & Stencil( int i = 0 ) {
    if(BaseGraph_.RowMap().GlobalIndicesInt())
      return RowStencil_int_[i];
    else
      throw "EpetraExt::BlockCrsMatrix::Stencil: Global Indices not int";
  }

  //! RowIndex
  int RowIndex( int i = 0 ) {
    if(BaseGraph_.RowMap().GlobalIndicesInt())
      return RowIndices_int_[i];
    else
      throw "EpetraExt::BlockCrsMatrix::RowIndex: Global Indices not int";
  }
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  //! Local Stencil Info
  const std::vector<long long> & Stencil64( int i = 0 ) {
    if(BaseGraph_.RowMap().GlobalIndicesLongLong())
      return RowStencil_LL_[i];
    else
      throw "EpetraExt::BlockCrsMatrix::Stencil: Global Indices not long long";
  }

  //! RowIndex
  long long RowIndex64( int i = 0 ) {
    if(BaseGraph_.RowMap().GlobalIndicesLongLong())
      return RowIndices_LL_[i];
    else
      throw "EpetraExt::BlockCrsMatrix::RowIndex: Global Indices not long long";
  }
#endif

  //! Routine for loading a base matrices values into the large Block Matrix
  //! The Row and Col arguments are indices into RowStencil 
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  void LoadBlock(const Epetra_RowMatrix & BaseMatrix, const int Row, const int Col);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  void LoadBlock(const Epetra_RowMatrix & BaseMatrix, const long long Row, const long long Col);
#endif

  //! Routine for summing base matrices values into the large Block Matrix
  //! The Row and Col arguments are indices into RowStencil 
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  void SumIntoBlock(double alpha, const Epetra_RowMatrix & BaseMatrix, const int Row, const int Col);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  void SumIntoBlock(double alpha, const Epetra_RowMatrix & BaseMatrix, const long long Row, const long long Col);
#endif

  //! Routine for summing base matrices values into the large Block Matrix
  //! The Row and Col arguments are global indices
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  void SumIntoGlobalBlock(double alpha, const Epetra_RowMatrix & BaseMatrix, const int Row, const int Col);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  void SumIntoGlobalBlock(double alpha, const Epetra_RowMatrix & BaseMatrix, const long long Row, const long long Col);
#endif

  //! Sum Entries into Block matrix using base-matrix numbering plus block Row and Col
  //! The Row and Col arguments are indices into RowStencil 
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  void BlockSumIntoGlobalValues(const int BaseRow, int NumIndices,
     double* Values, const int* Indices, const int Row, const int Col);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  void BlockSumIntoGlobalValues(const long long BaseRow, int NumIndices,
     double* Values, const long long* Indices, const long long Row, const long long Col);
#endif

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  void BlockReplaceGlobalValues(const int BaseRow, int NumIndices,
     double* Values, const int* Indices, const int Row, const int Col);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  void BlockReplaceGlobalValues(const long long BaseRow, int NumIndices,
     double* Values, const long long* Indices, const long long Row, const long long Col);
#endif

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  void BlockExtractGlobalRowView(const int BaseRow, int& NumEntries, 
     double*& Values, const int Row, const int Col);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  void BlockExtractGlobalRowView(const long long BaseRow, int& NumEntries, 
     double*& Values, const long long Row, const long long Col);
#endif

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  void ExtractBlock(Epetra_CrsMatrix & BaseMatrix, const int Row, const int Col);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  void ExtractBlock(Epetra_CrsMatrix & BaseMatrix, const long long Row, const long long Col);
#endif

 protected:

  Epetra_CrsGraph BaseGraph_;

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  std::vector< std::vector<int> > RowStencil_int_;

  std::vector<int> RowIndices_int_; 
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  std::vector< std::vector<long long> > RowStencil_LL_;

  std::vector<long long> RowIndices_LL_; 
#endif

  long long ROffset_;
  long long COffset_;

private:
  template<typename int_type>
  void TLoadBlock(const Epetra_RowMatrix & BaseMatrix, const int_type Row, const int_type Col);

  template<typename int_type>
  void TSumIntoBlock(double alpha, const Epetra_RowMatrix & BaseMatrix, const int_type Row, const int_type Col);

  template<typename int_type>
  void TSumIntoGlobalBlock(double alpha, const Epetra_RowMatrix & BaseMatrix, const int_type Row, const int_type Col);

  template<typename int_type>
  void TBlockSumIntoGlobalValues(const int_type BaseRow, int NumIndices,
     double* Values, const int_type* Indices, const int_type Row, const int_type Col);

  template<typename int_type>
  void TBlockReplaceGlobalValues(const int_type BaseRow, int NumIndices,
     double* Values, const int_type* Indices, const int_type Row, const int_type Col);

  template<typename int_type>
  void TBlockExtractGlobalRowView(const int_type BaseRow, int& NumEntries, 
     double*& Values, const int_type Row, const int_type Col);

  template<typename int_type>
  void TExtractBlock(Epetra_CrsMatrix & BaseMatrix, const int_type Row, const int_type Col);

  template<typename int_type>
  std::vector< std::vector<int_type> >& TRowStencil();

  template<typename int_type>
  std::vector<int_type>& TRowIndices();
};

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
template<> inline std::vector<int>& BlockCrsMatrix::TRowIndices() { return RowIndices_int_; }
template<> inline std::vector< std::vector<int> >& BlockCrsMatrix::TRowStencil() { return RowStencil_int_; }
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
template<> inline std::vector<long long>& BlockCrsMatrix::TRowIndices() { return RowIndices_LL_; }
template<> inline std::vector< std::vector<long long> >& BlockCrsMatrix::TRowStencil() { return RowStencil_LL_; }
#endif

} //namespace EpetraExt

#endif /* EPETRA_CRSMATRIX_H */
