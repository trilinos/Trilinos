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

//! EpetraExt::BlockUtility: A class of utilities for constructing block data structs

namespace EpetraExt {

class BlockUtility {
 public:

  //! BlockCrsMatrix constuctor with one block row per processor.
  /*! Creates a BlockGraph object
    
	\param In
	BaseGraph - Graph determining individual block structure, can be distrib. over subset of proc.'s
	\param In 
	RowStencil - Describes the stencil for block row on this processor (i.e. (-1 0 1) centered difference)
	\param In
	RowIndex - Defines the index used for this block row.
  */
  static Epetra_CrsGraph * GenerateBlockGraph( const Epetra_CrsGraph & BaseGraph, const std::vector<int> & RowStencil, int RowIndex, const Epetra_Comm & GlobalComm );
  
};

} //namespace EpetraExt

#endif /* EPETRA_CRSMATRIX_H */
