/*
// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
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
*/

#ifndef _EPETRAEXT_REDISTOR_H_
#define _EPETRAEXT_REDISTOR_H_

#include "Amesos_ConfigDefs.h"

#include "Epetra_DistObject.h"

//! EpetraExt_Redistor: simple class to facilitate redistribution of \c Epetra_DistObject's.
/*! EpetraExt_Redistor: a simple class to facilitate the redistribution
  of \c Epetra_DistObject's. An example of use of this class is the
  redistribution of matrix, LHS and RHS from a given map, to another
  map, so that the given DistObject can be converted to a non-epetra
  format:
  - objects on non-linear maps can be mapped to linear maps;
  - objects distributed on N-processors can be redistributed on a linear map on M-processors.
*/

class EpetraExt_Redistor
{
  
public:

  //! Initializes object from \c SourceMap, to map with all elements on processor 0.
  EpetraExt_Redistor(const Epetra_BlockMap & SourceMap) 
  {
    EpetraExt_Redistor(SourceMap,1);
  }
  
  //! Initializes object from \c SourceMap, to linear map defined on the first \c NumTargetProcs processors.
  /*! Initializes object for redistributions from \c SourceMap, to a linear map
      (with the same number of elements), defined on the first \c NumTargetProcs.
  */
  EpetraExt_Redistor(const Epetra_BlockMap & SourceMap, const int NumTargetProcs);

  //! Initializes object from \c SourceMap to \c TargetMap.
  EpetraExt_Redistor(const Epetra_BlockMap & SourceMap, const Epetra_BlockMap & TargetMap);
  
  ~EpetraExt_Redistor();

  //! Imports distributed object, defined on \c SourceMap, from an object defined on \c TargetMap.
  inline int SourceImport(Epetra_DistObject & Source, Epetra_DistObject & Target) const {
    return( SourceImport(Source,Target,Insert) );
  }

  //! Imports distributed object, defined on \c TargetMap, from an object defined on \c SourceMap.
  inline int TargetImport(Epetra_DistObject & Source, Epetra_DistObject & Target) const 
  {
    return( TargetImport(Source,Target,Insert) );
  }

  //! Imports distributed object, defined on \c SourceMap, from an object defined on \c TargetMap using specified \c CombineMode.
  int SourceImport(Epetra_DistObject & Source, Epetra_DistObject & Target,
		   Epetra_CombineMode CombineMode) const;

  //! Imports distributed object, defined on \c TargetMap, from an object defined on \c SourceMap using specified \c CombineMode.
  int TargetImport(Epetra_DistObject & Source, Epetra_DistObject & Target,
		   Epetra_CombineMode CombineMode) const;

  //! Returns a pointer to the \c SourceMap.
  inline Epetra_BlockMap * SourceMap() const
  {
    return( SourceMap_ );
  }

  //! Returns a pointer to the \c TargetMap.
  inline Epetra_BlockMap * TargetMap() const
  {
    return( TargetMap_ );
  }

  //! Returns a pointer to the \c Epetra_Import object from \c SourceMap to \c TargetMap.
  inline Epetra_Import * SourceImporter() const 
  {
    return( SourceImporter_ );
  }

  //! Returns a pointer to the \c Epetra_Import object from \c TargetMap to \c SourceMap.
  inline Epetra_Import * TargetImporter() const 
  {
    return( TargetImporter_ );
  }
  
private:

  Epetra_Import * SourceImporter_;       // import from SourceMap to TargetMap
  Epetra_Import * TargetImporter_;     

  Epetra_BlockMap * SourceMap_;
  Epetra_BlockMap * TargetMap_;
  
};

#endif
