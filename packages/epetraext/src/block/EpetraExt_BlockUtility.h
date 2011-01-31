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

#ifndef EPETRAEXT_BLOCKUTILITY_H
#define EPETRAEXT_BLOCKUTILITY_H

#include <vector>

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
  static Epetra_Map * GenerateBlockMap( const Epetra_BlockMap & BaseMap, const int*  RowIndices, int num_indices, const Epetra_Comm & GlobalComm );

  /*! Creates a BlockMap object
    
	\param In
	BaseMap - Map determining individual block structure, can be distrib. over subset of proc.'s
	\param In
	RowIndices - Defines the indices for local block rows
  */
  static Epetra_Map * GenerateBlockMap( const Epetra_BlockMap & BaseMap, const std::vector<int> & RowIndices, const Epetra_Comm & GlobalComm );

  /*! Creates a BlockMap object
    
	\param In
	BaseMap - Map determining individual block structure, can be distrib. over subset of proc.'s
	\param In
	RowIndices - Defines the indices for local block rows
  */
  static Epetra_Map * GenerateBlockMap( const Epetra_BlockMap & BaseMap, const Epetra_BlockMap& BlockMap, const Epetra_Comm & GlobalComm );

  //! BlockCrsMatrix constuctor
  /*! Creates a BlockGraph object
    
	\param In
	BaseGraph - Graph determining individual block structure, can be distrib. over subset of proc.'s
	\param In 
	RowStencil - Describes the stencil for block rows on this processor (i.e. (-1 0 1) centered difference)
	\param In
	RowIndices - Defines the indices for local block rows
  */
  static Epetra_CrsGraph * GenerateBlockGraph( const Epetra_CrsGraph & BaseGraph, const std::vector< std::vector<int> > & RowStencil, const std::vector<int> & RowIndices, const Epetra_Comm & GlobalComm );

  // Nearly identical version yet using RowMatrix interface instead of CrsGraph
  static Epetra_CrsGraph * GenerateBlockGraph( const Epetra_RowMatrix & BaseMatrix, const std::vector< std::vector<int> > & RowStencil, const std::vector<int> & RowIndices, const Epetra_Comm & GlobalComm );

  //! Generate global block graph using base graph and local block graph
  static Epetra_CrsGraph * GenerateBlockGraph( const Epetra_CrsGraph & BaseGraph, const Epetra_CrsGraph& LocalBlockGraph, const Epetra_Comm & GlobalComm );

  //! Generate stencil arrays from a local block graph
  static void GenerateRowStencil(const Epetra_CrsGraph& LocalBlockGraph, std::vector<int> RowIndices, std::vector< std::vector<int> >& RowStencil);

  //! Routine for calculating Offset for creating unique global IDs for Block representation
  static int CalculateOffset(const Epetra_BlockMap & BaseMap);
  
};

} //namespace EpetraExt

#endif /* EPETRA_CRSMATRIX_H */
