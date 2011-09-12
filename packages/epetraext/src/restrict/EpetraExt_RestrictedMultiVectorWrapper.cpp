/*
//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

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
