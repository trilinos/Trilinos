
#include <EpetraExt_Reindex_LinearProblem.h>

#include <EpetraExt_Reindex_CrsMatrix.h>
#include <EpetraExt_Reindex_MultiVector.h>

#include <Epetra_LinearProblem.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Map.h>

namespace EpetraExt {

LinearProblem_Reindex::
~LinearProblem_Reindex()
{
  if( newObj_ ) delete newObj_;

  if( MatTrans_ ) delete MatTrans_;
  if( LHSTrans_ ) delete LHSTrans_;
  if( RHSTrans_ ) delete RHSTrans_;

  if( NewRowMapOwned_ ) delete NewRowMap_;
}

LinearProblem_Reindex::NewTypeRef
LinearProblem_Reindex::
operator()( OriginalTypeRef orig )
{
  Epetra_CrsMatrix * OldMatrix = dynamic_cast<Epetra_CrsMatrix*>( orig.GetMatrix() );
  Epetra_MultiVector * OldRHS = orig.GetRHS();
  Epetra_MultiVector * OldLHS = orig.GetLHS();
  const Epetra_BlockMap & OldRowMap = OldMatrix->Map();

  //If new map not passed in create one with lex ordering
  if( !NewRowMap_ )
  {
    int NumMyElements = OldRowMap.NumMyElements();
    int NumGlobalElements = OldRowMap.NumGlobalElements();

    NewRowMap_ = new Epetra_Map( NumGlobalElements, NumMyElements, 0, OldRowMap.Comm() );
    NewRowMapOwned_ = true;
  }

  MatTrans_ = new CrsMatrix_Reindex( *NewRowMap_ );
  LHSTrans_ = new MultiVector_Reindex( *NewRowMap_ );
  RHSTrans_ = new MultiVector_Reindex( *NewRowMap_ );

  Epetra_CrsMatrix * NewMatrix = &((*MatTrans_)( *OldMatrix ));
  Epetra_MultiVector * NewLHS = &((*LHSTrans_)( *OldLHS ));
  Epetra_MultiVector * NewRHS = &((*RHSTrans_)( *OldRHS ));

  newObj_ = new Epetra_LinearProblem( NewMatrix, NewLHS, NewRHS );

  return *newObj_;
}

} //namespace EpetraExt

