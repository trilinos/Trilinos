
#include <EDT_MultiVector_View.h>

#include <Epetra_MultiVector.h>
#include <Epetra_BlockMap.h>

namespace Epetra_Transform {

std::auto_ptr<Epetra_MultiVector> MultiVector_View::operator()( const Epetra_MultiVector & original )
{
  int numVec = NumVec_;
  if( numVec == -1 ) numVec = original.NumVectors();

  double ** ptrArray;
  original.ExtractView( &ptrArray );
  return std::auto_ptr<Epetra_MultiVector>( new Epetra_MultiVector( View,
                                                                    NewMap_,
                                                                    ptrArray,
                                                                    numVec ) );
}

} //namespace Epetra_Transform
