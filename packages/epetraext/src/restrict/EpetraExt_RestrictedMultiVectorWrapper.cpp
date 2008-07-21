#include "EpetraExt_ConfigDefs.h"


#ifdef HAVE_MPI


#include "EpetraExt_RestrictedMultiVectorWrapper.h"
#include "Epetra_MpiComm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_MultiVector.h"


namespace EpetraExt {


RestrictedMultiVectorWrapper::RestrictedMultiVectorWrapper()
  : proc_is_active(true),
    subcomm_is_set(false),
    RestrictedComm_(0),
    MPI_SubComm_(MPI_COMM_NULL),
    ResMap_(0),
    input_mv_(),
    restricted_mv_()  
    {}


RestrictedMultiVectorWrapper::~RestrictedMultiVectorWrapper(){
  delete ResMap_;
  delete RestrictedComm_;
}


int RestrictedMultiVectorWrapper::SetMPISubComm(MPI_Comm MPI_SubComm){
  if(!subcomm_is_set){
    MPI_SubComm_=MPI_SubComm; delete RestrictedComm_; subcomm_is_set=true;
    return 0;
  }
  else return -1;
}


int RestrictedMultiVectorWrapper::restrict_comm(Teuchos::RCP<Epetra_MultiVector> input_mv){
 input_mv_=input_mv;

  /* Pull the Input Matrix Info */
  const Epetra_MpiComm *InComm = dynamic_cast<const Epetra_MpiComm*>(& input_mv_->Comm());
  const Epetra_BlockMap *InMap = dynamic_cast<const Epetra_BlockMap*>(& input_mv_->Map());

  if(!InComm || !InMap) return -1;
  
  if(!subcomm_is_set){
    /* Build the Split Communicators, If Needed */
    int color;
    if(InMap->NumMyElements()) color=1;
    else color=MPI_UNDEFINED;
    MPI_Comm_split(InComm->Comm(),color,InComm->MyPID(),&MPI_SubComm_);
  }
  else{
    /* Sanity check user-provided subcomm - drop an error if the MPISubComm
       does not include a processor with data. */
    if (input_mv->MyLength() && MPI_SubComm_ == MPI_COMM_NULL)
      return(-2);
  }

  /* Mark active processors */
  if(MPI_SubComm_ == MPI_COMM_NULL) proc_is_active=false;
  else proc_is_active=true;

  
  int Nrows=InMap->NumGlobalElements();
      
  if(proc_is_active){      
    RestrictedComm_=new Epetra_MpiComm(MPI_SubComm_);
    
    /* Build the Restricted Maps */
    ResMap_ = new Epetra_BlockMap(Nrows,InMap->NumMyElements(),InMap->MyGlobalElements(),
                                  InMap->ElementSizeList(),InMap->IndexBase(),*RestrictedComm_);

    /* Allocate the restricted matrix*/
    double *A; int LDA;
    input_mv_->ExtractView(&A,&LDA);
    restricted_mv_ = Teuchos::rcp(new Epetra_MultiVector(View,*ResMap_,A,LDA,input_mv_->NumVectors()));
  }
}/*end restrict_comm*/

}
  
#endif
