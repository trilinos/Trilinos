#include "EpetraExt_ConfigDefs.h"


#ifdef HAVE_MPI


#include "EpetraExt_RestrictedCrsMatrixWrapper.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"


namespace EpetraExt{

RestrictedCrsMatrixWrapper::RestrictedCrsMatrixWrapper()
  : proc_is_active(true),
    subcomm_is_set(false),
    MPI_SubComm_(MPI_COMM_NULL),
    RestrictedComm_(0),
    ResRowMap_(0),
    ResColMap_(0),
    input_matrix_(),
    restricted_matrix_()  
{
  
}

RestrictedCrsMatrixWrapper::~RestrictedCrsMatrixWrapper() {
  delete ResRowMap_;
  delete ResColMap_;
  delete RestrictedComm_;
}



int RestrictedCrsMatrixWrapper::SetMPISubComm(MPI_Comm MPI_SubComm){
  if(!subcomm_is_set){
    MPI_SubComm_=MPI_SubComm; delete RestrictedComm_; subcomm_is_set=true;
    return 0;
  }
  else return -1;
}



int RestrictedCrsMatrixWrapper::restrict_comm(Teuchos::RCP<Epetra_CrsMatrix> input_matrix){
  /* Pull the Matrix Info */
  input_matrix_=input_matrix;
  
  const Epetra_MpiComm *InComm = dynamic_cast<const Epetra_MpiComm*>(& input_matrix_->Comm());
  const Epetra_Map *InRowMap= dynamic_cast<const Epetra_Map* >(& input_matrix_->RowMap());
  const Epetra_Map *InColMap= dynamic_cast<const Epetra_Map* >(& input_matrix_->ColMap());

  if(!InComm || !InRowMap || !InColMap) return (-1);
  
  int Nrows=InRowMap->NumGlobalElements();
  int Ncols=InColMap->NumGlobalElements();
  
  if(!subcomm_is_set){
    /* Build the Split Communicators, If Needed */
    int color;
    if(InRowMap->NumMyElements()) color=1;
    else color=MPI_UNDEFINED;
    MPI_Comm_split(InComm->Comm(),color,InComm->MyPID(),&MPI_SubComm_);
  }
  else{
    /* Sanity check user-provided subcomm - drop an error if the MPISubComm
       does not include a processor with data. */
    if (input_matrix->NumMyRows() && MPI_SubComm_ == MPI_COMM_NULL)
      return(-2);
  }

  /* Mark active processors */
  if(MPI_SubComm_ == MPI_COMM_NULL) proc_is_active=false;
  else proc_is_active=true;
  

  if(proc_is_active){      
    RestrictedComm_=new Epetra_MpiComm(MPI_SubComm_);
    
    /* Build the Restricted Maps */
    ResRowMap_ = new Epetra_Map(Nrows,InRowMap->NumMyElements(),InRowMap->MyGlobalElements(),
                                InRowMap->IndexBase(),*RestrictedComm_);
    ResColMap_ = new Epetra_Map(Ncols,InColMap->NumMyElements(),InColMap->MyGlobalElements(),
                                InColMap->IndexBase(),*RestrictedComm_);        
    
    int *rowptr,*colind,Nr;
    double *values;
    
    /* Allocate the Restricted Matrix */
    restricted_matrix_= Teuchos::rcp(new Epetra_CrsMatrix(View,*ResRowMap_,*ResColMap_,0));
    for(int i=0;i<input_matrix_->NumMyRows();i++) {
      input_matrix_->ExtractMyRowView(i,Nr,values,colind);
      restricted_matrix_->InsertMyValues(i,Nr,values,colind);
    }
    restricted_matrix_->FillComplete();      
  }

  return 0;
}/*end restrict_comm*/



}  
#endif
