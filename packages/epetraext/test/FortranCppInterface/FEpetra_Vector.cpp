
#include "FEpetra_Vector.h"
#include "FEpetra_Vector_Cpp.hpp"
#include "FEpetra_Map_Cpp.hpp"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"


namespace {


using Teuchos::Array;
using Teuchos::RCP;


// Note that by storing RCP objects that when the program goes down, any
// vectors not destoryed explicitly will be destroyed here so valgrind will
// not complain!  I don't know if this is a good thing or a bad thing?

Array<RCP<Epetra_Vector> >& tableOfVectors()
{
  static Array<RCP<Epetra_Vector> > loc_tableOfVectors;
  return loc_tableOfVectors;
}
// 2007/12/18: rabartl: ToDo: Validate that the above table is empty when the
// program finishes!  


} // namespace


//
// Definitions from FEpetra_Vector.h
//


extern "C" {


VectorID FEpetra_Vector_Create( MapID mapID )
{
  using Teuchos::rcp;
  const RCP<const Epetra_Map> map = FEpetra::getMap(mapID);
  tableOfVectors().push_back(rcp(new Epetra_Vector(*map)));
  return tableOfVectors().size() - 1;
}


void FEpetra_Vector_Destroy( VectorID vectorID )
{
  tableOfVectors()[vectorID] = Teuchos::null;
  // 2007/12/18: rabartl: ToDo: Maintain a free list of vectorIDs that can be
  // reused.  This will allow for many many allocations and deallocations!
}


void FEpetra_Vector_Random( VectorID vectorID )
{
  FEpetra::getVector(vectorID)->Random();
}


void FEpetra_Vector_Update(
  VectorID vectorID, double alpha, VectorID vector2ID, double beta
  )
{
  FEpetra::getVector(vectorID)->Update(
    alpha, *FEpetra::getVector(vector2ID), beta );
}


double FEpetra_Vector_Norm2( VectorID vectorID )
{
  double norm2;
  FEpetra::getVector(vectorID)->Norm2(&norm2);
  return norm2;
}


} // extern "C"


//
// Definitions from FEpetra_Vector_Cpp.hpp
//


const Teuchos::RCP<Epetra_Vector>
FEpetra::getVector( VectorID vectorID )
{
  return tableOfVectors()[vectorID];
}
