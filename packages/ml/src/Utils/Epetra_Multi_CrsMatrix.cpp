#include "ml_config.h"
#include "Epetra_Multi_CrsMatrix.h"
#include "ml_epetra.h"
#include "ml_epetra_utils.h"
#include "ml_mat_formats.h"
#if defined(HAVE_ML_EPETRA)

#include "Epetra_Comm.h"
//#include "EpetraExt_MatrixMatrix.h"//haq

extern void Epetra_CrsMatrix_Print(const Epetra_CrsMatrix& A, ostream& os);//haq


void ML_Matrix_Print(ML_Operator *ML,const Epetra_Comm &Comm,const Epetra_Map &Map, char *fname){
  Epetra_CrsMatrix *Temp_;
  ofstream ofs(fname);
  Epetra_CrsMatrix_Wrap_ML_Operator(ML,Comm,Map,&Temp_,View);
  Epetra_CrsMatrix_Print(*Temp_,ofs);
  delete Temp_;
}







// ================================================ ====== ==== ==== == = 
// Constructor
ML_Epetra::Epetra_Multi_CrsMatrix::Epetra_Multi_CrsMatrix(int NumMatrices,Epetra_CrsMatrix ** CrsMatrices)
  :NumMatrices_(NumMatrices),CrsMatrices_(CrsMatrices),Label_(0)
{
  Label_=strdup("Epetra_Multi_CrsMatrix");
  Comm_ = &(CrsMatrices_[0]->Comm());
  DomainMap_ = &(CrsMatrices_[NumMatrices_-1]->OperatorDomainMap());
  RangeMap_ = &(CrsMatrices_[0]->OperatorRangeMap());
}/*end Constructor*/

// ================================================ ====== ==== ==== == = 
// Destructor
ML_Epetra::Epetra_Multi_CrsMatrix::~Epetra_Multi_CrsMatrix(){if(Label_) free(Label_);}


// ================================================ ====== ==== ==== == = 
// Computes Y= <me> * X
int ML_Epetra::Epetra_Multi_CrsMatrix::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  int nv=X.NumVectors();
  Epetra_MultiVector *MV[2]={0,0};
  MV[(NumMatrices_-1)%2]=(Epetra_MultiVector*)&X;

  /* Do the first n-1 matvecs */
  for(int i=NumMatrices_-1;i>0;i--){
    //    printf("[%d] MCRS: matvec %d begins\n",Comm_->MyPID(),i);
    if(MV[(i+1)%2] && MV[(i+1)%2]!=&X) delete MV[(i+1)%2];
    MV[(i+1)%2]=new Epetra_MultiVector(CrsMatrices_[i]->RangeMap(),nv,false);
    //    printf("[%d] MCRS: Matvec %d / %d ov=%d iv=%d mat=%dx%d\n",Comm_->MyPID(),i,NumMatrices_,           
    //           MV[(i+1)%2]->GlobalLength(),MV[i%2]->GlobalLength(),
    //           CrsMatrices_[i]->NumGlobalRows(),CrsMatrices_[i]->NumGlobalCols());
    CrsMatrices_[i]->Apply(*MV[i%2],*MV[(i+1)%2]);
    //    printf("[%d] MCRS: matvec complete\n",Comm_->MyPID());
  }/*end for*/

  /* Final matvec */
  //  printf("[%d] MCRS: Matvec %d / %d \n",Comm_->MyPID(),0,NumMatrices_);  
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
  //  printf("[%d] Epetra_Multi_CrsMatrix::MatrixMatrix_Multiply[ML] called\n",A.Comm().MyPID());
  ML_Comm* temp = global_comm;
  
  /* Setup for 1st Matmat */
  ML_Operator * MV[2]={0,0},*CV;
  MV[(NumMatrices_-1)%2]= ML_Operator_Create(comm);  
  //  printf("[%d] MATMAT: Prewrap called\n",A.Comm().MyPID());
  rv=ML_Operator_WrapEpetraCrsMatrix((Epetra_CrsMatrix*)&A,MV[(NumMatrices_-1)%2]);
  ML_CHK_ERR(rv);

  /* Do the matmats */
  for(int i=NumMatrices_-1;rv==0 && i>=0;i--){
    //    printf("[%d] Matmat %d/%d Starting\n",A.Comm().MyPID(),NumMatrices_-i,NumMatrices_);
    if(MV[(i+1)%2] && i!=NumMatrices_-1) ML_Operator_Destroy(&MV[(i+1)%2]);
    MV[(i+1)%2]=ML_Operator_Create(comm);
    CV=ML_Operator_Create(comm);
    rv=ML_Operator_WrapEpetraCrsMatrix(CrsMatrices_[i],CV);
    ML_CHK_ERR(rv);


    //    printf("[%d] Matmat %d/%d Prewrapped\n",A.Comm().MyPID(),NumMatrices_-i,NumMatrices_);
    //    printf("[%d] MV[i%%2]=%#x MV[(i+1)%%2]=%#x\n",A.Comm().MyPID(),MV[i%2],MV[(i+1)%2]);

    ML_2matmult(CV,MV[i%2],MV[(i+1)%2],ML_CSR_MATRIX);

    /* DEBUG */
    char str[80];
    sprintf(str,"op11.%d.dat",NumMatrices_-1-i);
    ML_Matrix_Print(MV[(i+1)%2],A.Comm(),CrsMatrices_[i]->RangeMap(),str);




    ML_Operator_Destroy(&CV);
    //    printf("[%d] Matmat %d/%d Complete\n",A.Comm().MyPID(),NumMatrices_-i,NumMatrices_);    
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
  //  printf("[%d] Epetra_Multi_CrsMatrix::MatrixMatrix_Multiply[Epetra] called\n",A.Comm().MyPID());
  ML_Comm* comm;
  ML_Comm_Create(&comm);
  ML_Operator *C_;
  int rv=MatrixMatrix_Multiply(A,comm,&C_);
  Epetra_CrsMatrix_Wrap_ML_Operator(C_,*Comm_,*RangeMap_,C,Copy);
  ML_Operator_Destroy(&C_);
  ML_Comm_Destroy(&comm);
  return rv;
}/*end MatrixMatrix_Multiply*/
#endif
