// @HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef EpetraExt_CRSGRAPH_SUBCOPY_H
#define EpetraExt_CRSGRAPH_SUBCOPY_H

#include <EpetraExt_Transform.h>
#include <Epetra_Map.h>

class Epetra_CrsMatrix;
class Epetra_Import;

namespace EpetraExt {

//! Generates a sub-block view of a Epetra_CrsMatrix
class CrsMatrix_SubCopy : public SameTypeTransform<Epetra_CrsMatrix> {

  Epetra_Map  newRowMap_;
  Epetra_Map  newColMap_;
  Epetra_Map  newDomainMap_;
  Epetra_Map  newRangeMap_;
  Epetra_Import * importer_;


 public:

  //! Destructor
  ~CrsMatrix_SubCopy();

  //! Constructor
  CrsMatrix_SubCopy( const Epetra_Map  & newMap )
  : newRowMap_(newMap),
    newColMap_(newMap),
    newDomainMap_(newMap),
    newRangeMap_(newMap)
  {}

  //! Constructor
  CrsMatrix_SubCopy( const Epetra_Map  & newRangeAndRowMap,
		    const Epetra_Map  & newDomainMap )
  : newRowMap_(newRangeAndRowMap),
    newColMap_(newDomainMap),
    newDomainMap_(newDomainMap),
    newRangeMap_(newRangeAndRowMap)
  {}

  //! Transformation Operator
  NewTypeRef operator()( OriginalTypeRef orig );

  ///
  bool fwd();

  ///
  bool rvs();
};

} //namespace EpetraExt

#endif //EpetraExt_CRSGRAPH_SUBCOPY_H
