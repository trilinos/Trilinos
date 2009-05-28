#include "Ifpack_HIPS.h"
#include "Ifpack_Utils.h"
extern "C" {
#include "hips.h"
};
//#include "io.h"
#include "Epetra_MpiComm.h"
#include "Epetra_IntVector.h"
#include "Epetra_Import.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#ifdef IFPACK_NODE_AWARE_CODE
extern int ML_NODE_ID;
#endif

using Teuchos::RefCountPtr;
using Teuchos::rcp;

#if defined(HAVE_IFPACK_HIPS) && defined(HAVE_MPI)


Ifpack_HIPS::Ifpack_HIPS(Epetra_RowMatrix* A):
  A_(rcp(A,false)),
  IsParallel_(false),
  IsInitialized_(false),
  IsComputed_(false),
  Label_(),
  HIPS_id(-1), // Assumes user will initialze HIPS outside 
  Time_(A_->Comm())
{
  
}

void Ifpack_HIPS::Destroy(){
  //NTS: Assume user will call HIPS_Finalize elsewhere - HIPS_Clean never needs
  //to be called if HIPS_Finalize is called at the end, unless you want to reuse
  //a slot.
}



int Ifpack_HIPS::Initialize(){
  if(Comm().NumProc() != 1) IsParallel_ = true;
  else IsParallel_ = false;
  IsInitialized_=true;
  return 0;
}

int Ifpack_HIPS::SetParameters(Teuchos::ParameterList& parameterlist){
  List_=parameterlist;

  // Grab the hips ID
  HIPS_id=List_.get("hips: id",-1);
  if(HIPS_id==-1) IFPACK_CHK_ERR(-1);

  // Set Defaults
  HIPS_SetDefaultOptions(HIPS_id,List_.get("hips: strategy",HIPS_ITERATIVE));  
  
  // Set the communicator
  const Epetra_MpiComm* MpiComm=dynamic_cast<const Epetra_MpiComm*>(&A_->Comm());
  if(!MpiComm) IFPACK_CHK_ERR(-2);  
  HIPS_SetCommunicator(HIPS_id,MpiComm->GetMpiComm());
  
  // Options
  HIPS_SetOptionINT(HIPS_id,HIPS_SYMMETRIC,List_.get("hips: symmetric",0));
  HIPS_SetOptionINT(HIPS_id,HIPS_VERBOSE,List_.get("hips: setup output",1));
  HIPS_SetOptionINT(HIPS_id,HIPS_LOCALLY,List_.get("hips: fill",0));
  // Sadly, this fill option doesn't work for HIPS_ITERATIVE mode, meaning the
  // only way to control fill-in is via the drop tolerance. 
  HIPS_SetOptionINT(HIPS_id,HIPS_ITMAX,List_.get("hips: krylov its",1));  
  HIPS_SetOptionINT(HIPS_id,HIPS_KRYLOV_METHOD,List_.get("hips: krylov",0));
  HIPS_SetOptionINT(HIPS_id,HIPS_REORDER,List_.get("hips: reorder",1));
  HIPS_SetOptionINT(HIPS_id,HIPS_GRAPH_SYM,List_.get("hips: graph symmetric",0));
  HIPS_SetOptionINT(HIPS_id,HIPS_FORTRAN_NUMBERING,List_.get("hips: fortran numbering",0));
  HIPS_SetOptionINT(HIPS_id,HIPS_DOF,List_.get("hips: dof per node",1));

  // Make sure the ILU always runs, by setting the internal tolerance really, really low.
  HIPS_SetOptionREAL(HIPS_id,HIPS_PREC,1e-100);

  // Other HAX
  HIPS_SetOptionINT(HIPS_id,HIPS_KRYLOV_RESTART,1);    
  
  // Options for Iterative only
  if(List_.get("hips: strategy",HIPS_ITERATIVE)==HIPS_ITERATIVE){  
    HIPS_SetOptionREAL(HIPS_id, HIPS_DROPTOL0, List_.get("hips: drop tolerance",1e-2));
    HIPS_SetOptionREAL(HIPS_id, HIPS_DROPTOL1, List_.get("hips: drop tolerance",1e-2));
    HIPS_SetOptionREAL(HIPS_id, HIPS_DROPTOLE, List_.get("hips: drop tolerance",1e-2));    
  }
  // NTS: This is only a subset of the actual HIPS options. 
  return 0;
}


int Ifpack_HIPS::Compute(){
  if(HIPS_id==-1) IFPACK_CHK_ERR(-1);
  int N=A_->NumMyRows(), nnz=A_->NumMyNonzeros();

  // Pull the column indices, if possible
  int *rowptr,*colind,ierr,maxnr,Nr;
  double *values;  
  Epetra_CrsMatrix *Acrs=dynamic_cast<Epetra_CrsMatrix*>(&*A_);
  const Epetra_Map &RowMap=A_->RowMatrixRowMap();
  const Epetra_Map &ColMap=A_->RowMatrixColMap();
  if(Acrs) Acrs->ExtractCrsDataPointers(rowptr,colind,values);
  else{
    maxnr=A_->MaxNumEntries();
    colind=new int[maxnr];
    values=new double[maxnr];
  }

  // Calculate MinGIDs (HIPS likes things starting at 0/1
  // NTS: This is necessary for the node filter stuff
  int MinGRID=RowMap.MinAllGID();
  int MinGCID=ColMap.MinAllGID();
  
  // Build a new row map w/ MinGID 0 (for HIPS)
  const int * OldGIDs=RowMap.MyGlobalElements();  
  int *newGIDs=new int[N];
  for(int i=0;i<N;i++)
    newGIDs[i]=OldGIDs[i]-MinGRID; 
  RowMap0_=rcp(new Epetra_Map(-1,N,newGIDs,0,A_->Comm()));
  delete [] newGIDs;

  
  //NTS: Need to totally rewrite the graph insertion algorithms since they're
  //gratuitous when we already have a matrix assembled. I'd do
  //it now, but I want "functional" first and "not terribly slow" second
  ierr=HIPS_GraphBegin(HIPS_id,A_->NumGlobalRows(),nnz);
  if(ierr!=HIPS_SUCCESS) IFPACK_CHK_ERR(-2);  
  if(Acrs){
    // Graph insert - CRS mode
    for(int i=0;i<N;i++)
      for(int j=rowptr[i];j<rowptr[i+1];j++){
        ierr=HIPS_GraphEdge(HIPS_id,Acrs->GRID(i)-MinGRID,Acrs->GCID(colind[j])-MinGCID);
        if(ierr!=HIPS_SUCCESS) IFPACK_CHK_ERR(-3);
      }
  }
  else{
    // Graph insert - RM mode
    for(int i=0;i<N;i++){
      A_->ExtractMyRowCopy(i,maxnr,Nr,values,colind);
      for(int j=0;j<Nr;j++){
        ierr=HIPS_GraphEdge(HIPS_id,RowMap.GID(i)-MinGRID,ColMap.GID(colind[j])-MinGCID);
        if(ierr!=HIPS_SUCCESS) IFPACK_CHK_ERR(-3);
      }
    }    
  }
  ierr=HIPS_GraphEnd(HIPS_id);
  if(ierr!=HIPS_SUCCESS) IFPACK_CHK_ERR(-4);  
  
  /*Have processor 0 send in the partition*/
  // NTS: This is really, really annoying.  Look at all this import/export
  // stuff.  This is mind-numbingly unnecessary.
  const Epetra_Comm &Comm=A_->Comm();
  int mypid=Comm.MyPID();
  
  Epetra_Map OnePerProcMap(-1,1,0,Comm);
  Epetra_IntVector RowsPerProc(OnePerProcMap);
  Epetra_IntVector RowGID(View,*RowMap0_,RowMap0_->MyGlobalElements());
  
  // Get the RPP partial sums
  Comm.ScanSum(&N,&(RowsPerProc[0]),1); 

  // Build the maps for xfer to proc 0
  int OPP_els=0,RPP_els=0;
  if(!mypid){OPP_els=Comm.NumProc(); RPP_els=A_->NumGlobalRows();}
  Epetra_Map OPPMap_0(-1,OPP_els,0,Comm);
  Epetra_Map RPPMap_0(-1,RPP_els,0,Comm);
  Epetra_Import OPP_importer(OPPMap_0,OnePerProcMap);
  Epetra_Import RPP_importer(RPPMap_0,*RowMap0_);
  
  // Pull the vectors to proc 0
  Epetra_IntVector OPP_0(OPPMap_0);
  Epetra_IntVector RPP_0(RPPMap_0);
  OPP_0.Import(RowsPerProc,OPP_importer,Add);
  RPP_0.Import(RowGID,RPP_importer,Add);
 
  // Setup the partition
  if(!mypid){    
    int *mapptr=0;
    mapptr=new int[Comm.NumProc()+1];
    mapptr[0]=0;
    for(int i=0;i<Comm.NumProc();i++)
      mapptr[i+1]=OPP_0[i];

    // Call is only necessary on proc 0
    ierr=HIPS_SetPartition(HIPS_id,A_->Comm().NumProc(),mapptr,RPP_0.Values());
    HIPS_ExitOnError(ierr);
    delete [] mapptr;
  }      
  
  ierr = HIPS_AssemblyBegin(HIPS_id, nnz, HIPS_ASSEMBLY_OVW, HIPS_ASSEMBLY_OVW, HIPS_ASSEMBLY_FOOL,0);
  if(ierr!=HIPS_SUCCESS){HIPS_PrintError(ierr);IFPACK_CHK_ERR(-5);}

  if(Acrs){
    // Matrix insert - CRS Mode
    for(int i=0;i<N;i++){     
      for(int j=rowptr[i];j<rowptr[i+1];j++){      
        ierr = HIPS_AssemblySetValue(HIPS_id,Acrs->GRID(i)-MinGRID,Acrs->GCID(colind[j])-MinGCID, values[j]);
        if(ierr!=HIPS_SUCCESS){HIPS_PrintError(ierr);IFPACK_CHK_ERR(-6);}
      }
    }
  }
  else{
    // Matrix insert - RM Mode
    for(int i=0;i<N;i++){
      A_->ExtractMyRowCopy(i,maxnr,Nr,values,colind);
      for(int j=0;j<Nr;j++){      
        ierr = HIPS_AssemblySetValue(HIPS_id,RowMap.GID(i)-MinGRID,ColMap.GID(colind[j])-MinGCID, values[j]);
        if(ierr!=HIPS_SUCCESS){HIPS_PrintError(ierr);IFPACK_CHK_ERR(-6);}
      }
    }
  }  
  ierr = HIPS_AssemblyEnd(HIPS_id);
  if(ierr!=HIPS_SUCCESS){HIPS_PrintError(ierr);IFPACK_CHK_ERR(-7);}

  // Force factorization
  //NTS: This is odd.  There should be a better way to force this to happen.
  //  Epetra_Vector X(A_->RowMap()),Y(A_->RowMap());
  //  X.PutScalar(1.0);
  double *X=new double[3*N];
  double *Y=new double[3*N];
  for(int i=0;i<3*N;i++) X[i]=1.0;

  // NTS: Fix this to follow the ApplyInverse stuff so we're guaranteed not to
  // step on memory.
  
  ierr=HIPS_SetLocalRHS(HIPS_id,&(X[0]),HIPS_ASSEMBLY_OVW,HIPS_ASSEMBLY_OVW);
  if(ierr!=HIPS_SUCCESS) IFPACK_CHK_ERR(-11);

  ierr=HIPS_GetLocalSolution(HIPS_id,&(Y[0]));
  if(ierr!=HIPS_SUCCESS) IFPACK_CHK_ERR(-12);

  // Reset output for iteration
  HIPS_SetOptionINT(HIPS_id,HIPS_VERBOSE,List_.get("hips: iteration output",0));

  // Set Label
  int nnzP=0;
  HIPS_GetInfoINT(HIPS_id,HIPS_INFO_NNZ,&nnzP);
  if(nnzP>0) sprintf(Label_,"Ifpack_HIPS [dt=%4.1e fill=%4.2f]",List_.get("hips: drop tolerance",1e-2),(double)nnzP/(double)A_->NumGlobalNonzeros());
  else sprintf(Label_,"Ifpack_HIPS [dt=%4.1e]",List_.get("hips: drop tolerance",1e-2));
  // NTS: fill requires a HIPS debug level of at least 2
  
  IsComputed_=true;

  // Cleanup
  if(!Acrs){
    delete [] colind;
    delete [] values;
  }
  delete [] X;
  delete [] Y;
  return 0;
}




int Ifpack_HIPS::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const{
  int rv;
  if (!IsComputed())
    IFPACK_CHK_ERR(-3);

  if (X.NumVectors() != Y.NumVectors())
    IFPACK_CHK_ERR(-2);

  // HAQ: For now
  if(X.NumVectors()!=1) IFPACK_CHK_ERR(-42);

  // Wrapping for X==Y
  Teuchos::RefCountPtr< Epetra_MultiVector > X2;
  if (X.Pointers()[0] == Y.Pointers()[0])
    X2 = Teuchos::rcp( new Epetra_MultiVector(X) );    
  else
    X2 = Teuchos::rcp( (Epetra_MultiVector*)&X, false );

  Time_.ResetStartTime();

  // Force HIPS to do it's own Import/Export
  rv=HIPS_SetRHS(HIPS_id,A_->NumMyRows(),RowMap0_->MyGlobalElements(),(*X2)[0],HIPS_ASSEMBLY_OVW,HIPS_ASSEMBLY_OVW,HIPS_ASSEMBLY_FOOL);
  if(rv!=HIPS_SUCCESS) IFPACK_CHK_ERR(-11);
  
  rv=HIPS_GetSolution(HIPS_id,A_->NumMyRows(),RowMap0_->MyGlobalElements(),Y[0],HIPS_ASSEMBLY_FOOL);
  if(rv!=HIPS_SUCCESS) {
    HIPS_PrintError(rv);
    IFPACK_CHK_ERR(-12);
  }
  
  ++NumApplyInverse_;
  ApplyInverseTime_ += Time_.ElapsedTime();

  return(0);
}


ostream& Ifpack_HIPS::Print(ostream& os) const{
  os<<"Need to add meaningful output"<<endl;
  return os;
}


double Ifpack_HIPS::Condest(const Ifpack_CondestType CT, 
                             const int MaxIters,
                             const double Tol,
                             Epetra_RowMatrix* Matrix_in){
  return -1.0;
}




#endif
