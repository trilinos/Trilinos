
#include <EpetraExt_SolverMap_LinearProblem.h>

#include <Epetra_LinearProblem.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>

class Epetra_MultiVector;

namespace EpetraExt {

LinearProblem_SolverMap::
~LinearProblem_SolverMap()
{
  if( newObj_ && newObj_ != origObj_ ) delete newObj_;
}

LinearProblem_SolverMap::NewTypeRef
LinearProblem_SolverMap::
operator()( OriginalTypeRef orig )
{
  origObj_ = &orig;

  Epetra_CrsMatrix * OldMatrix = dynamic_cast<Epetra_CrsMatrix*>( orig.GetMatrix() );
  Epetra_MultiVector * OldRHS = orig.GetRHS();
  Epetra_MultiVector * OldLHS = orig.GetLHS();

  Epetra_CrsMatrix & NewMatrix = SolverMapTrans_( *OldMatrix );

  if( &NewMatrix == OldMatrix ) //same matrix so use same problem
    newObj_ = origObj_;
  else
    newObj_ = new Epetra_LinearProblem( &NewMatrix, OldLHS, OldRHS );

  return *newObj_;
}

} //namespace EpetraExt

