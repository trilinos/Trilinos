
#include <EDT_MultiVector_View.h>

#include <Epetra_MultiVector.h>
#include <Epetra_BlockMap.h>

namespace EpetraExt {

MultiVector_View::
~MultiVector_View()
{
  if( newObj_ ) delete newObj_;
}

MultiVector_View::NewTypeRef
MultiVector_View::
operator()( OriginalTypeRef orig )
{
  origObj_ = &orig;

  int numVec = NumVec_;
  if( numVec == -1 ) numVec = orig.NumVectors();

  double ** ptrArray;
  orig.ExtractView( &ptrArray );

  Epetra_MultiVector * newMV = new Epetra_MultiVector( View, NewMap_, ptrArray, numVec );

  newObj_ = newMV;

  return *newMV;
}

} // namespace EpetraExt

