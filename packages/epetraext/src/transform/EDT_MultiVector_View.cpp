
#include <EDT_MultiVector_View.h>

#include <Epetra_MultiVector.h>
#include <Epetra_BlockMap.h>

EpetraExt::MultiVector_View::NewTypePtr EpetraExt::MultiVector_View::operator()( EpetraExt::MultiVector_View::OriginalTypeRef original )
{
  int numVec = NumVec_;
  if( numVec == -1 ) numVec = original.NumVectors();

  double ** ptrArray;
  original.ExtractView( &ptrArray );
  return new Epetra_MultiVector( View, NewMap_, ptrArray, numVec );
}

