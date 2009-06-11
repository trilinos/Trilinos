/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */
#include "ml_config.h"
#if defined(HAVE_ML_EPETRA)
#include "Epetra_Multi_CrsMatrix.h"
#include "ml_epetra.h"
#include "ml_epetra_utils.h"
#include "ml_mat_formats.h"

#include "Epetra_Comm.h"
#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#endif

// ================================================ ====== ==== ==== == = 
// Constructor
ML_Epetra::Epetra_Multi_CrsMatrix::Epetra_Multi_CrsMatrix(int NumMatrices,Epetra_CrsMatrix ** CrsMatrices)
  :NumMatrices_(NumMatrices),CrsMatrices_(CrsMatrices),Label_(0)
{
  Label_=new char[80];
  strcpy(Label_,"Epetra_Multi_CrsMatrix");
  Comm_ = &(CrsMatrices_[0]->Comm());
  DomainMap_ = &(CrsMatrices_[NumMatrices_-1]->OperatorDomainMap());
  RangeMap_ = &(CrsMatrices_[0]->OperatorRangeMap());
}/*end Constructor*/

// ================================================ ====== ==== ==== == = 
// Destructor
ML_Epetra::Epetra_Multi_CrsMatrix::~Epetra_Multi_CrsMatrix(){if(Label_) delete [] Label_;}


// ================================================ ====== ==== ==== == = 
// Computes Y= <me> * X
int ML_Epetra::Epetra_Multi_CrsMatrix::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  int nv=X.NumVectors();
  Epetra_MultiVector *MV[2]={0,0};
  MV[(NumMatrices_-1)%2]=(Epetra_MultiVector*)&X;

  /* Do the first n-1 matvecs */
  for(int i=NumMatrices_-1;i>0;i--){
    if(MV[(i+1)%2] && MV[(i+1)%2]!=&X) delete MV[(i+1)%2];
    MV[(i+1)%2]=new Epetra_MultiVector(CrsMatrices_[i]->RangeMap(),nv,false);
    CrsMatrices_[i]->Apply(*MV[i%2],*MV[(i+1)%2]);
  }/*end for*/

  /* Final matvec */
  CrsMatrices_[0]->Apply(*MV[0],Y);    

  /* Cleanup */
  if(MV[1] != &X) delete MV[1];
  if(MV[0] != &X) delete MV[0];
  return 0;
}/*end Apply*/

// ================================================ ====== ==== ==== == = 
// Computes C= <me> * A
int  ML_Epetra::Epetra_Multi_CrsMatrix::MatrixMatrix_Multiply(const Epetra_CrsMatrix & A, ML_Comm *comm, ML_Operator **C) const
{

  int rv=0;
  ML_Comm* temp = global_comm;  
  /* DEBUG*/
  char str[80];

  /* Setup for 1st Matmat */
  ML_Operator * MV[2]={0,0},*CV;
  MV[(NumMatrices_-1)%2]= ML_Operator_Create(comm);
  rv=ML_Operator_WrapEpetraCrsMatrix((Epetra_CrsMatrix*)&A,MV[(NumMatrices_-1)%2]);
  ML_CHK_ERR(rv);

  /* Do the matmats */
  for(int i=NumMatrices_-1;rv==0 && i>=0;i--){
    /* Do pre-wraps */
    if(MV[(i+1)%2] && i!=NumMatrices_-1) ML_Operator_Destroy(&MV[(i+1)%2]);
    MV[(i+1)%2]=ML_Operator_Create(comm);
    CV=ML_Operator_Create(comm);
    rv=ML_Operator_WrapEpetraCrsMatrix(CrsMatrices_[i],CV);
    ML_CHK_ERR(rv);

    /* Do matmat */
    ML_2matmult(CV,MV[i%2],MV[(i+1)%2],ML_CSR_MATRIX);

    ML_Operator_Destroy(&CV);
  }/*end for*/
  global_comm = temp;

  /* Final Postwrap */
  *C=MV[1];

  /* Cleanup */
  if(MV[0]) ML_Operator_Destroy(&MV[0]);

  return rv;
}/*end MatrixMatrix_Multiply*/


// ================================================ ====== ==== ==== == = 
// Computes C= <me> * A
int ML_Epetra::Epetra_Multi_CrsMatrix::MatrixMatrix_Multiply(const Epetra_CrsMatrix & A, Epetra_CrsMatrix **C) const
{
  ML_Comm* comm;
  ML_Comm_Create(&comm);
#ifdef ML_MPI
  const Epetra_MpiComm *epcomm = dynamic_cast<const Epetra_MpiComm*>(&(A.Comm()));
  // Get the MPI communicator, as it may not be MPI_COMM_W0RLD, and update the ML comm object
  if (epcomm) ML_Comm_Set_UsrComm(comm,epcomm->Comm());
#endif
  ML_Operator *C_;
  int rv=MatrixMatrix_Multiply(A,comm,&C_);
  Epetra_CrsMatrix_Wrap_ML_Operator(C_,*Comm_,*RangeMap_,C,Copy);//,A.IndexBase());
  ML_Operator_Destroy(&C_);
  ML_Comm_Destroy(&comm);
  return rv;
}/*end MatrixMatrix_Multiply*/
#endif
