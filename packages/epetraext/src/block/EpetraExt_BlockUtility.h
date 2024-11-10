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

#ifndef EPETRAEXT_BLOCKUTILITY_H
#define EPETRAEXT_BLOCKUTILITY_H

#if defined(EpetraExt_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The EpetraExt package is deprecated"
#endif
#endif

#include <vector>

#include <Epetra_ConfigDefs.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_RowMatrix.h>

//! EpetraExt::BlockUtility: A class of utilities for constructing block data structs

namespace EpetraExt {

class BlockUtility {
 public:

  /*! Creates a BlockMap object
    
	\param In
	BaseMap - Map determining individual block structure, can be distrib. over subset of proc.'s
	\param In
	RowIndices - Defines the indices for local block rows
  */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  static Epetra_Map * GenerateBlockMap( const Epetra_BlockMap & BaseMap, const int*  RowIndices, int num_indices, const Epetra_Comm & GlobalComm, int Offset = 0 );
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  static Epetra_Map * GenerateBlockMap( const Epetra_BlockMap & BaseMap, const long long*  RowIndices, int num_indices, const Epetra_Comm & GlobalComm, long long Offset = 0 );
#endif

  /*! Creates a BlockMap object
    
	\param In
	BaseMap - Map determining individual block structure, can be distrib. over subset of proc.'s
	\param In
	RowIndices - Defines the indices for local block rows
  */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  static Epetra_Map * GenerateBlockMap( const Epetra_BlockMap & BaseMap, const std::vector<int> & RowIndices, const Epetra_Comm & GlobalComm, int Offset = 0 );
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  static Epetra_Map * GenerateBlockMap( const Epetra_BlockMap & BaseMap, const std::vector<long long> & RowIndices, const Epetra_Comm & GlobalComm, long long Offset = 0 );
#endif

  /*! Creates a BlockMap object
    
	\param In
	BaseMap - Map determining individual block structure, can be distrib. over subset of proc.'s
	\param In
	RowIndices - Defines the indices for local block rows
  */
  static Epetra_Map * GenerateBlockMap( const Epetra_BlockMap & BaseMap, const Epetra_BlockMap& BlockMap, const Epetra_Comm & GlobalComm, int Offset = 0 );

  //! BlockCrsMatrix constuctor
  /*! Creates a BlockGraph object
    
	\param In
	BaseGraph - Graph determining individual block structure, can be distrib. over subset of proc.'s
	\param In 
	RowStencil - Describes the stencil for block rows on this processor (i.e. (-1 0 1) centered difference)
	\param In
	RowIndices - Defines the indices for local block rows
  */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  static Epetra_CrsGraph * GenerateBlockGraph( const Epetra_CrsGraph & BaseGraph, const std::vector< std::vector<int> > & RowStencil, const std::vector<int> & RowIndices, const Epetra_Comm & GlobalComm );
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  static Epetra_CrsGraph * GenerateBlockGraph( const Epetra_CrsGraph & BaseGraph, const std::vector< std::vector<long long> > & RowStencil, const std::vector<long long> & RowIndices, const Epetra_Comm & GlobalComm );
#endif

  // Nearly identical version yet using RowMatrix interface instead of CrsGraph
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  static Epetra_CrsGraph * GenerateBlockGraph( const Epetra_RowMatrix & BaseMatrix, const std::vector< std::vector<int> > & RowStencil, const std::vector<int> & RowIndices, const Epetra_Comm & GlobalComm );
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  static Epetra_CrsGraph * GenerateBlockGraph( const Epetra_RowMatrix & BaseMatrix, const std::vector< std::vector<long long> > & RowStencil, const std::vector<long long> & RowIndices, const Epetra_Comm & GlobalComm );
#endif

  //! Generate global block graph using base graph and local block graph
  static Epetra_CrsGraph * GenerateBlockGraph( const Epetra_CrsGraph & BaseGraph, const Epetra_CrsGraph& LocalBlockGraph, const Epetra_Comm & GlobalComm );

  //! Generate stencil arrays from a local block graph
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  static void GenerateRowStencil(const Epetra_CrsGraph& LocalBlockGraph, std::vector<int> RowIndices, std::vector< std::vector<int> >& RowStencil);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  static void GenerateRowStencil(const Epetra_CrsGraph& LocalBlockGraph, std::vector<long long> RowIndices, std::vector< std::vector<long long> >& RowStencil);
#endif

  //! Routine for calculating Offset for creating unique global IDs for Block representation
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  static int CalculateOffset(const Epetra_BlockMap & BaseMap);
#endif

  static long long CalculateOffset64(const Epetra_BlockMap & BaseMap);

private:

  template<typename int_type>
  static int_type TCalculateOffset(const Epetra_BlockMap & BaseMap);

  template<typename int_type>
  static Epetra_Map * TGenerateBlockMap( const Epetra_BlockMap & BaseMap, const int_type*  RowIndices, int num_indices, const Epetra_Comm & GlobalComm, int_type Offset = 0);

  template<typename int_type>
  static Epetra_CrsGraph * TGenerateBlockGraph( const Epetra_CrsGraph & BaseGraph, const std::vector< std::vector<int_type> > & RowStencil, const std::vector<int_type> & RowIndices, const Epetra_Comm & GlobalComm );

  template<typename int_type>
  static Epetra_CrsGraph * TGenerateBlockGraph( const Epetra_RowMatrix & BaseMatrix, const std::vector< std::vector<int_type> > & RowStencil, const std::vector<int_type> & RowIndices, const Epetra_Comm & GlobalComm );

  template<typename int_type>
  static Epetra_CrsGraph * TGenerateBlockGraph( const Epetra_CrsGraph & BaseGraph, const Epetra_CrsGraph& LocalBlockGraph, const Epetra_Comm & GlobalComm );

  template<typename int_type>
  static void TGenerateRowStencil(const Epetra_CrsGraph& LocalBlockGraph, std::vector<int_type> RowIndices, std::vector< std::vector<int_type> >& RowStencil);
};

} //namespace EpetraExt

#endif /* EPETRA_CRSMATRIX_H */
