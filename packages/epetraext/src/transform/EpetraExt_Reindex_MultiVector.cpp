
#include <EpetraExt_Reindex_MultiVector.h>

#include <vector>

#include <Epetra_Map.h>
#include <Epetra_MultiVector.h>

namespace EpetraExt {

MultiVector_Reindex::
~MultiVector_Reindex()
{
  if( newObj_ ) delete newObj_;
}

MultiVector_Reindex::NewTypeRef
MultiVector_Reindex::
operator()( OriginalTypeRef orig )
{
  origObj_ = &orig;

  //test map, must have same number of local and global elements as original row map
  assert( orig.Map().NumMyElements() == NewRowMap_.NumMyElements() );

  vector<double*> MyValues(1);
  int MyLDA;
  int NumVectors = orig.NumVectors();
  orig.ExtractView( &MyValues[0], &MyLDA );

  Epetra_MultiVector * NewMV = new Epetra_MultiVector( View, NewRowMap_, MyValues[0], MyLDA, NumVectors );

  newObj_ = NewMV;

  return *NewMV;
}

} // namespace EpetraExt

