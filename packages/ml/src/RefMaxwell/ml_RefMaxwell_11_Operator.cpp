/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */
#include "ml_config.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_EPETRAEXT)
#include "ml_epetra.h"
#include "ml_epetra.h"
#include "ml_RefMaxwell_11_Operator.h"
#include "ml_epetra_utils.h"
#include "EpetraExt_Transpose_RowMatrix.h"
#include "EpetraExt_SolverMap_CrsMatrix.h"
#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#endif


#define NO_OUTPUT

/* Disabling this option should result in a slight speedup, but will break matvecs with the 11-Operator */
#define MANUALLY_TRANSPOSE_D0

#ifdef NO_OUTPUT
#define ML_Matrix_Print(w,x,y,z) ;
#else
extern void ML_Matrix_Print(ML_Operator *ML,const Epetra_Comm &Comm,const Epetra_Map &Map, char *fname);
#endif

// ================================================ ====== ==== ==== == = 
// Constructor
ML_Epetra::ML_RefMaxwell_11_Operator::ML_RefMaxwell_11_Operator(const Epetra_CrsMatrix& SM_Matrix,    //S+M
                                                             const Epetra_CrsMatrix& D0_Matrix,    //T or D0
                                                             const Epetra_CrsMatrix& M0inv_Matrix, //M0^{-1}
                                                             const Epetra_CrsMatrix& M1_Matrix):   //M1(1)
  SM_Matrix_(&SM_Matrix),Addon_Matrix_(0),D0T_Matrix_(0)
{
  Label_=new char [80];
  strcpy(Label_,"ML_RefMaxwell_11_Operator");
  Comm_ = &(SM_Matrix_->Comm());
  DomainMap_ = &(SM_Matrix_->OperatorDomainMap());
  RangeMap_ = &(SM_Matrix_->OperatorRangeMap());

  /* Transpose D0 */
#ifdef MANUALLY_TRANSPOSE_D0
  D0_Matrix_Transposer_= new EpetraExt::RowMatrix_Transpose((Epetra_Map*)&M0inv_Matrix.OperatorRangeMap());
  D0T_Matrix_= dynamic_cast<Epetra_CrsMatrix*>( & ((*D0_Matrix_Transposer_)((Epetra_CrsMatrix&)D0_Matrix)));
  D0T_Matrix_= dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(*D0T_Matrix_,D0T_Matrix_Trans_,"D0T",false));
#endif
  
  /* Build the Epetra_Multi_CrsMatrix */
  Addon_Matrix_=new Epetra_CrsMatrix*[5];
  Addon_Matrix_[0]=Addon_Matrix_[4]=(Epetra_CrsMatrix*)&M1_Matrix;
  Addon_Matrix_[1]=(Epetra_CrsMatrix*)&D0_Matrix;
#ifdef MANUALLY_TRANSPOSE_D0
  Addon_Matrix_[3]=D0T_Matrix_;
#else
  Addon_Matrix_[3]=(Epetra_CrsMatrix*)&D0_Matrix;//HAQ
#endif
  Addon_Matrix_[2]=(Epetra_CrsMatrix*)&M0inv_Matrix;
  Addon_=new Epetra_Multi_CrsMatrix(5,Addon_Matrix_);

}/*end Constructor*/ 

// ================================================ ====== ==== ==== == = 
// Destructor
ML_Epetra::ML_RefMaxwell_11_Operator::~ML_RefMaxwell_11_Operator()
{
  if(Label_) delete [] Label_;
  if(Addon_Matrix_) delete [] Addon_Matrix_;
  if(Addon_) delete Addon_;
#ifdef MANUALLY_TRANSPOSE_D0
  if(D0_Matrix_Transposer_) delete D0_Matrix_Transposer_;
#endif
}/*end destructor*/


// ================================================ ====== ==== ==== == = 
// Computes Y= <me> * X
int ML_Epetra::ML_RefMaxwell_11_Operator::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  Epetra_MultiVector temp(X);
  /* Do the SM part */
  SM_Matrix_->Apply(X,Y);

  /* Do the Addon part */
  Addon_->Apply(X,temp);

  /* Sum things together*/
  Y.Update(1,temp,1);
  return 0;
}/*end Apply*/


// ================================================ ====== ==== ==== == = 
// Computes C= A^T * <me> * A.  OptimizeStorage *must* be called for both A and the
// matrices in *this, before this routine can work.
int ML_Epetra::ML_RefMaxwell_11_Operator::PtAP(const Epetra_CrsMatrix & P, ML_Comm *comm, ML_Operator **C) const{
  ML_Operator *SM_ML,*P_ML,*R_ML,*PtSMP_ML,*temp1,*temp2,*opwrap,*D0_M1_P_ML;

  /* General Stuff */  
  ML_Comm* temp = global_comm;
  P_ML  = ML_Operator_Create(comm);
  R_ML  = ML_Operator_Create(comm);;
  ML_Operator_WrapEpetraCrsMatrix((Epetra_CrsMatrix*)&P,P_ML);  
  ML_Operator_Transpose_byrow(P_ML,R_ML);

  /* Do the SM part */
  SM_ML = ML_Operator_Create(comm);
  temp1 = ML_Operator_Create(comm);
  PtSMP_ML  = ML_Operator_Create(comm);
  ML_Operator_WrapEpetraCrsMatrix((Epetra_CrsMatrix*)SM_Matrix_,SM_ML);
  ML_2matmult(SM_ML,P_ML,temp1,ML_CSR_MATRIX);
  ML_2matmult_block(R_ML,temp1,PtSMP_ML,ML_CSR_MATRIX);
  ML_Operator_Destroy(&temp1);
  ML_Operator_Destroy(&SM_ML);
#ifdef MANUALLY_TRANSPOSE_D0
  ML_Operator_Destroy(&R_ML);
#endif
  ML_Matrix_Print(PtSMP_ML,*Comm_,*RangeMap_,"ptsmp.dat");

#ifdef MANUALLY_TRANSPOSE_D0
  /* Do the Addon: Step #1: M1 * P*/
  opwrap = ML_Operator_Create(comm);
  temp1 = ML_Operator_Create(comm);
  ML_Operator_WrapEpetraCrsMatrix((Epetra_CrsMatrix*)Addon_Matrix_[4],opwrap);
  ML_2matmult(opwrap,P_ML,temp1,ML_CSR_MATRIX);
  ML_Operator_Destroy(&opwrap);
  ML_Operator_Destroy(&P_ML);

  /* Do the Addon: Step #2: D0^T *(M1 * P)*/
  opwrap = ML_Operator_Create(comm);
  D0_M1_P_ML = ML_Operator_Create(comm);
  ML_Operator_WrapEpetraCrsMatrix((Epetra_CrsMatrix*)Addon_Matrix_[3],opwrap);
  ML_2matmult(opwrap,temp1,D0_M1_P_ML,ML_CSR_MATRIX);
  ML_Operator_Destroy(&opwrap);
  ML_Operator_Destroy(&temp1);
  
  /* Do the Addon: Step #3: M0^{-1} * (D0^T * M1 * P)*/
  opwrap = ML_Operator_Create(comm);
  temp1 = ML_Operator_Create(comm);
  ML_Operator_WrapEpetraCrsMatrix((Epetra_CrsMatrix*)Addon_Matrix_[2],opwrap);
  ML_2matmult(opwrap,D0_M1_P_ML,temp1,ML_CSR_MATRIX);
  ML_Operator_Destroy(&opwrap);  

  /* Do the Addon: Step #4: Transpose (D0^T * M1 * P) & multiply by output from Step 3*/
  opwrap = ML_Operator_Create(comm);
  temp2 = ML_Operator_Create(comm);  
  ML_Operator_Transpose_byrow(D0_M1_P_ML,opwrap);
  ML_2matmult(opwrap,temp1,temp2,ML_CSR_MATRIX);
  ML_Operator_Destroy(&opwrap);
  ML_Operator_Destroy(&temp1);
  ML_Operator_Destroy(&D0_M1_P_ML);
  ML_Matrix_Print(temp2,*Comm_,*RangeMap_,"pt_add_p.dat");
#else
  ML_Operator *P_M1_D0_ML;

  /* Do the Addon: Step #1: P^T * M1 */
  opwrap = ML_Operator_Create(comm);
  temp1 = ML_Operator_Create(comm);
  ML_Operator_WrapEpetraCrsMatrix((Epetra_CrsMatrix*)Addon_Matrix_[0],opwrap);
  ML_2matmult_block(R_ML,opwrap,temp1,ML_CSR_MATRIX);
  ML_Operator_Destroy(&opwrap);
  ML_Operator_Destroy(&P_ML);
  ML_Operator_Destroy(&R_ML);
  
  /* Do the Addon: Step #2: (P^T * M1) * D0*/
  opwrap = ML_Operator_Create(comm);
  P_M1_D0_ML = ML_Operator_Create(comm);
  ML_Operator_WrapEpetraCrsMatrix((Epetra_CrsMatrix*)Addon_Matrix_[1],opwrap);
  ML_2matmult_block(temp1,opwrap,P_M1_D0_ML,ML_CSR_MATRIX);
  ML_Operator_Destroy(&opwrap);
  ML_Operator_Destroy(&temp1);  

  /* Do the Addon: Step #3: (P^T * M1 * D0) * M0^{-1} */
  opwrap = ML_Operator_Create(comm);
  temp1 = ML_Operator_Create(comm);
  ML_Operator_WrapEpetraCrsMatrix((Epetra_CrsMatrix*)Addon_Matrix_[2],opwrap);
  ML_2matmult(P_M1_D0_ML,opwrap,temp1,ML_CSR_MATRIX);
  ML_Operator_Destroy(&opwrap);  

  /* Do the Addon: Step #4: Transpose (P^T * M1 * D0) & multiply by output from Step 3*/
  opwrap = ML_Operator_Create(comm);
  temp2 = ML_Operator_Create(comm);  
  ML_Operator_Transpose_byrow(P_M1_D0_ML,opwrap);
  ML_2matmult(temp1,opwrap,temp2,ML_CSR_MATRIX);
  ML_Operator_Destroy(&opwrap);
  ML_Operator_Destroy(&temp1);
  ML_Operator_Destroy(&P_M1_D0_ML);
  ML_Matrix_Print(temp2,*Comm_,*RangeMap_,"pt_add_p_rev.dat");   
#endif

  /* Add the matrices together */
  ML_Operator_Add(PtSMP_ML,temp2,*C,ML_CSR_MATRIX,1.0);
  ML_Matrix_Print(*C,*Comm_,*RangeMap_,"ptap.dat");  

  /* Cleanup */
  global_comm = temp;
  ML_Operator_Destroy(&temp2);
  ML_Operator_Destroy(&PtSMP_ML);

  return 0;
}


// ================================================ ====== ==== ==== == = 
// Computes C= <me> * A
int  ML_Epetra::ML_RefMaxwell_11_Operator::MatrixMatrix_Multiply(const Epetra_CrsMatrix & A, ML_Comm *comm, ML_Operator **C) const
{
  ML_Operator *SM_ML,*A_ML,*temp1,*temp2;

  /* General Stuff */  
  ML_Comm* temp = global_comm;
  A_ML  = ML_Operator_Create(comm);
  *C = ML_Operator_Create(comm);
  ML_Operator_WrapEpetraCrsMatrix((Epetra_CrsMatrix*)&A,A_ML);  
  
  /* Do the SM part */
  SM_ML = ML_Operator_Create(comm);
  temp1 = ML_Operator_Create(comm);
  ML_Operator_WrapEpetraCrsMatrix((Epetra_CrsMatrix*)SM_Matrix_,SM_ML);
  ML_2matmult(SM_ML,A_ML,temp1,ML_CSR_MATRIX);
  ML_Matrix_Print(temp1,*Comm_,*RangeMap_,"smp.dat");

  /* Do the Addon part */
  Addon_->MatrixMatrix_Multiply(A,comm,&temp2);
  ML_Matrix_Print(temp2,*Comm_,*RangeMap_,"add_p.dat");
  
  /* Add the matrices together */
  ML_Operator_Add(temp2,temp1,*C,ML_CSR_MATRIX,1.0);
  ML_Matrix_Print(*C,*Comm_,*RangeMap_,"tfinal.dat");  
  
  /* Cleanup */
  global_comm = temp;
  ML_Operator_Destroy(&A_ML);
  ML_Operator_Destroy(&SM_ML);
  ML_Operator_Destroy(&temp1);
  ML_Operator_Destroy(&temp2);

  return 0;
}/*end MatrixMatrix_Multiply*/


// ================================================ ====== ==== ==== == = 
// Computes C= <me> * A
int ML_Epetra::ML_RefMaxwell_11_Operator::MatrixMatrix_Multiply(const Epetra_CrsMatrix & A, Epetra_CrsMatrix **C) const
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
  Epetra_CrsMatrix_Wrap_ML_Operator(C_,*Comm_,*DomainMap_,C,Copy,A.IndexBase());
  ML_Operator_Destroy(&C_);
  ML_Comm_Destroy(&comm);
  return rv;
}/*end MatrixMatrix_Multiply*/
#endif
