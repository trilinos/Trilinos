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

#include "Amesos_Klu.h"
#include "Amesos_Support.h"

Amesos_StandardIndex::Amesos_StandardIndex( const Epetra_Map& OriginalMap ) 
{

#ifdef HAVE_AMESOS_EPETRAEXT
  int NumGlobalElements = OriginalMap.NumGlobalElements();
  int NumMyElements = OriginalMap.NumMyElements();
  StdIndexMap_ = rcp( new Epetra_Map( NumGlobalElements, NumMyElements, 0, OriginalMap.Comm() ) );
  MatTrans_ = rcp( new EpetraExt::CrsMatrix_Reindex( *StdIndexMap_ ) );
  VecTrans_ = rcp( new EpetraExt::MultiVector_Reindex( *StdIndexMap_ ) );
#endif

}

#ifdef HAVE_AMESOS_EPETRAEXT
  //! Convert MultiVector to a MultiVector indexed from 0 to n-1 
Epetra_MultiVector* Amesos_StandardIndex::StandardizeIndex( Epetra_MultiVector* OriginalMultiVector ) {

  return ( &((*VecTrans_)( *OriginalMultiVector )) );
  
}

  //! Convert MultiVector to a MultiVector indexed from 0 to n-1 
Teuchos::RCP<Epetra_MultiVector> Amesos_StandardIndex::StandardizeIndex( Epetra_MultiVector & OriginalMultiVector ) {

  return VecTrans_->transform(OriginalMultiVector);
}
  

  //! Convert CrsMatrix to a CrsMatrix indexed from 0 to n-1 
Epetra_CrsMatrix*  Amesos_StandardIndex::StandardizeIndex( Epetra_CrsMatrix* OriginalCrsMatrix ) {

    return &((*MatTrans_)( *OriginalCrsMatrix ));
}
#endif

