#ifndef AMESOS_SUPPORT_H
#define AMESOS_SUPPORT_H

#if defined(Amesos_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Amesos package is deprecated"
#endif
#endif

#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"

#ifdef HAVE_AMESOS_EPETRAEXT
#include "EpetraExt_Reindex_CrsMatrix.h"
#include "EpetraExt_Reindex_MultiVector.h"
#endif
/*!
 \class Amesos_Support
 
 \brief Amesos_Support: Collection of utilities not included in Amesos.h

 \author Ken Stanley

 \date Last updated on 16-July-05
*/

#if 0
class Amesos_Support
{
public:
  //! Default constructor.
  Amesos_Support() {}
  
  //! Default destructor.
  ~Amesos_Support() {}

};
#endif

#if 0
//! Returns a Matrix indexed from 0 to n-1 
Epetra_CrsMatrix* Amesos_StandardIndexMatrix( const Epetra_CrsMatrix&* OriginalMatrix );
#endif

class Amesos_StandardIndex
{
 public:
  //! Default constructor.
  Amesos_StandardIndex( const Epetra_Map& OriginalMap ) ;
  
  //! Default destructor.
  ~Amesos_StandardIndex() {}

#ifdef HAVE_AMESOS_EPETRAEXT
  //! Convert MultiVector to a MultiVector indexed from 0 to n-1 
  Epetra_MultiVector* StandardizeIndex( Epetra_MultiVector* OriginalMultiVector );

  //! Convert MultiVector to a MultiVector indexed from 0 to n-1 
  Teuchos::RCP<Epetra_MultiVector> StandardizeIndex( Epetra_MultiVector & OriginalMultiVector );

  //! Convert CrsMatrix to a CrsMatrix indexed from 0 to n-1 
  Epetra_CrsMatrix* StandardizeIndex( Epetra_CrsMatrix* OriginalCrsMatrix );

  //! Convert CrsMatrix to a CrsMatrix indexed from 0 to n-1 
  Epetra_Map*  StdIndexMap() {
    return &*StdIndexMap_ ; 
  }
#endif
    


private:
#ifdef HAVE_AMESOS_EPETRAEXT
  //! Points to a Map which standardized indices - i.e. from 0 to n-1 
  Teuchos::RCP<Epetra_Map> StdIndexMap_;
  //! Points to an object which reindexes a CrsMatrix to a contiguous map
  Teuchos::RCP<EpetraExt::CrsMatrix_Reindex> MatTrans_;
  //! Points to an object which reindexes a MultiVector to a contiguous map
  Teuchos::RCP<EpetraExt::MultiVector_Reindex> VecTrans_;
#endif
} ;
#endif
