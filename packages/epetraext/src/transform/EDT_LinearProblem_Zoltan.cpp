
#include <EDT_LinearProblem_Zoltan.h>
#include <EDT_CrsGraph_Zoltan.h>

#include <Epetra_Export.h>
#include <Epetra_Import.h>
#include <Epetra_LinearProblem.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Epetra_IntVector.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>

namespace EpetraExt {

LinearProblem_Zoltan::~LinearProblem_Zoltan()
{
  if( Exporter_ ) delete Exporter_;
  if( Importer_ ) delete Importer_;

  if( NewProblem_ ) delete NewProblem_;
  if( NewRHS_ ) delete NewRHS_;
  if( NewLHS_ ) delete NewLHS_;
  if( NewMatrix_ ) delete NewMatrix_;
  if( NewGraph_ ) delete NewGraph_;
  if( NewRowMap_ ) delete NewRowMap_;

}

LinearProblem_Zoltan::NewTypePtr LinearProblem_Zoltan::operator()( LinearProblem_Zoltan::OriginalTypeRef original )
{
  OldProblem_ = &original;
  OldMatrix_ = dynamic_cast<Epetra_CrsMatrix*>( original.GetMatrix() );
  OldGraph_ = const_cast<Epetra_CrsGraph*>(&OldMatrix_->Graph());
  OldRHS_ = original.GetRHS();
  OldLHS_ = original.GetLHS();
  OldRowMap_ = const_cast<Epetra_Map*>(&OldMatrix_->RowMap());

  int ierr = 0;

  const Epetra_Comm & CommObj = OldRowMap_->Comm();

  if( !OldMatrix_ ) ierr = -2;
  if( !OldRHS_ )    ierr = -3;
  if( !OldLHS_ )    ierr = -4;

  NewGraph_ = CrsGraph_Zoltan()(*OldGraph_);
  NewMatrix_ = new Epetra_CrsMatrix( Copy, *NewGraph_ );

  Epetra_BlockMap * tmpMap = new Epetra_BlockMap( NewGraph_->RowMap() );
  NewRowMap_ = dynamic_cast<Epetra_Map*>(tmpMap);
  if( !NewRowMap_ ) abort();

  NewRHS_ = new Epetra_MultiVector( *NewRowMap_, 1 );
  NewLHS_ = new Epetra_MultiVector( *NewRowMap_, 1 );

  Exporter_ = new Epetra_Export( *OldRowMap_, *NewRowMap_ );
  Importer_ = new Epetra_Import( *NewRowMap_, *OldRowMap_ );

  NewProblem_ = new Epetra_LinearProblem( NewMatrix_, NewLHS_, NewRHS_ );

  return NewProblem_;
}

bool EpetraExt::LinearProblem_Zoltan::fwd()
{
  NewLHS_->Export( *OldLHS_, *Exporter_, Insert );
  NewRHS_->Export( *OldRHS_, *Exporter_, Insert );
  NewMatrix_->Export( *OldMatrix_, *Exporter_, Insert );

  return true;
}

bool EpetraExt::LinearProblem_Zoltan::rvs()
{
  OldLHS_->Import( *NewLHS_, *Exporter_, Insert );
  OldRHS_->Import( *NewRHS_, *Exporter_, Insert );
  OldMatrix_->Import( *NewMatrix_, *Exporter_, Insert );

  return true;
}

} //namespace EpetraExt
