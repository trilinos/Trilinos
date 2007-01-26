#include "ml_config.h"
#include "ml_RefMaxwell_11_Operator.h"
//#include "ml_epetra.h"
#include "ml_epetra_utils.h"
//#include "ml_mat_formats.h"
#include "EpetraExt_Transpose_RowMatrix.h"
#include "EpetraExt_SolverMap_CrsMatrix.h"

#include "EpetraExt_MatrixMatrix.h" //haq


#include "EpetraExt_RowMatrixOut.h"// debugging

#include "EpetraExt_Transpose_RowMatrix.h" //haq?
#include "EpetraExt_SolverMap_CrsMatrix.h" //haq?


#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_EPETRAEXT)


extern Epetra_RowMatrix* ModifyEpetraMatrixColMap(const Epetra_RowMatrix &A,
                                                  EpetraExt::CrsMatrix_SolverMap &transform);//haq

extern void print_stats(const Epetra_CrsMatrix& A, char *label);//haq
// ================================================ ====== ==== ==== == = 
// Constructor
ML_Epetra::ML_RefMaxwell_11_Operator::ML_RefMaxwell_11_Operator(const Epetra_CrsMatrix& SM_Matrix,    //S+M
                                                             const Epetra_CrsMatrix& D0_Matrix,    //T or D0
                                                             const Epetra_CrsMatrix& M0inv_Matrix, //M0^{-1}
                                                             const Epetra_CrsMatrix& M1_Matrix):   //M1(1)
  SM_Matrix_(&SM_Matrix),Addon_Matrix_(0),Core_Matrix_(0)
{
  Label_=strdup("ML_RefMaxwell_11_Operator");
  Comm_ = &(SM_Matrix_->Comm());
  DomainMap_ = &(SM_Matrix_->OperatorDomainMap());
  RangeMap_ = &(SM_Matrix_->OperatorRangeMap());
  //  printf("[%d] ML_RefMaxwell_11_Operator: Constructor Called\n",Comm_->MyPID());


  /* Transpose D0 */
  //  EpetraExt::RowMatrix_Transpose transposer(true,(Epetra_Map*)&M0inv_Matrix.OperatorRangeMap());
  //  Epetra_CrsMatrix &D0_Matrix_transpose= dynamic_cast<Epetra_CrsMatrix&>(transposer((Epetra_CrsMatrix&)D0_Matrix));

  /* Build the DO M0 D0^T Core */
  //  Epetra_CrsMatrix *core_temp[2] = {(Epetra_CrsMatrix*)&D0_Matrix,(Epetra_CrsMatrix*)&M0inv_Matrix};
  //  Epetra_Multi_CrsMatrix CoreTemp(2,core_temp);
  //  int rv=CoreTemp.MatrixMatrix_Multiply(D0_Matrix_transpose,&Core_Matrix_);a
  //  if(rv) printf("[%d] ERROR: Core Matrix formation failed with error code %d\n",Comm_->MyPID(),rv);
  //  else  printf("[%d] ML_RefMaxwell_11_Operator: Core Built\n",Comm_->MyPID());
  
  /* Build the D0 M0 D0^T Core with EpetraExt */
  Epetra_CrsMatrix temp1(Copy,M0inv_Matrix.RangeMap(),0);
  Core_Matrix_=new Epetra_CrsMatrix(Copy,*RangeMap_,0);
  
  EpetraExt::MatrixMatrix::Multiply(M0inv_Matrix,false,D0_Matrix,true,temp1);
  EpetraExt::MatrixMatrix::Multiply(D0_Matrix,false,temp1,false,*Core_Matrix_);

  //NTS: For whatever reason ML's wrap routines die *horribly* when we play with
  //the transposed matrix created by EpetraExt::RowMatrix_Transpose.  This
  //should be fixed eventually.  As for now, we'll just use EpetraExt to build
  //the core.
  
  //  printf("[%d] ML_RefMaxwell_11_Operator: MCRS Constructor Done\n",Comm_->MyPID());

  Core_Matrix_->OptimizeStorage();
  
  
  /* Build the Epetra_Multi_CrsMatrix */
  Addon_Matrix_=new Epetra_CrsMatrix*[3];
  Addon_Matrix_[0]=Addon_Matrix_[2]=(Epetra_CrsMatrix*)&M1_Matrix;
  Addon_Matrix_[1]=Core_Matrix_;
  Addon_=new Epetra_Multi_CrsMatrix(3,Addon_Matrix_);

  //  printf("[%d] ML_RefMaxwell_11_Operator: Multi_CrsMatrix Built\n",Comm_->MyPID());
  
}/*end Constructor*/ 

// ================================================ ====== ==== ==== == = 
// Destructor
ML_Epetra::ML_RefMaxwell_11_Operator::~ML_RefMaxwell_11_Operator()
{
  if(Label_) free(Label_);
  if(Addon_Matrix_) delete Addon_Matrix_;
  if(Addon_) delete Addon_;
  if(Core_Matrix_) delete Core_Matrix_;
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
// Computes C= <me> * A
int  ML_Epetra::ML_RefMaxwell_11_Operator::MatrixMatrix_Multiply(const Epetra_CrsMatrix & A, ML_Comm *comm, ML_Operator **C) const
{  
  ML_Operator *SM_ML,*A_ML,*temp1,*temp2;
  //  printf("[%d] RM11:MMM [ML]\n",A.Comm().MyPID());
  
  /* General Stuff */
  ML_Comm* temp = global_comm;
  A_ML  = ML_Operator_Create(comm);
  *C = ML_Operator_Create(comm);
  ML_Operator_WrapEpetraCrsMatrix((Epetra_CrsMatrix*)&A,A_ML);  
  
  /* Do the SM part */
  //  printf("[%d] RM11: SM Part of Matmat commencing\n",A.Comm().MyPID());
  SM_ML = ML_Operator_Create(comm);
  temp1 = ML_Operator_Create(comm);
  ML_Operator_WrapEpetraCrsMatrix((Epetra_CrsMatrix*)SM_Matrix_,SM_ML);
  ML_2matmult(SM_ML,A_ML,temp1,ML_CSR_MATRIX);


#define SANITY_CHECK 1
#ifdef SANITY_CHECK
  /* DEBUG */
  Epetra_CrsMatrix C_EP(Copy,*DomainMap_,0);
  EpetraExt::MatrixMatrix::Multiply(*SM_Matrix_,false,A,false,C_EP);
  Epetra_CrsMatrix diff(C_EP), *C_EP2;
  Epetra_CrsMatrix_Wrap_ML_Operator(temp1,*Comm_,*DomainMap_,&C_EP2,Copy);
  C_EP2->OptimizeStorage();

  
  //  printf("[%d] RM11: SM sanity check commencing\n",A.Comm().MyPID());
  EpetraExt::MatrixMatrix::Add(*C_EP2,false,1.0,diff,-1.0);

  
  double norm1=C_EP.NormInf();  
  double norm2=C_EP2->NormInf();
  double norm3=diff.NormInf()/norm1;
  if(A.Comm().MyPID()==0)
     printf("RM11: diff = %6.4e (%6.4e/%6.4e)\n",norm3,norm1,norm2);
  delete C_EP2;
  /*end DEBUG*/
#endif
  
  /* Do the Addon part */
  //  printf("[%d] RM11: Addon Part of Matmat commencing\n",A.Comm().MyPID());
  Addon_->MatrixMatrix_Multiply(A,comm,&temp2);

  //  ML_Operator_Print_UsingGlobalOrdering(temp2,"addon.ml",0,0);

  
  /* Add the matrices together */
  //  printf("[%d] RM11: Final Matrix Add commencing\n",A.Comm().MyPID());
  ML_Operator_Add(temp1,temp2,*C,ML_CSR_MATRIX,1.0);

  //  ML_Operator_Print_UsingGlobalOrdering(*C,"res.ml",0,0);
  
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
  //  printf("[%d] RM11:MMM [Epetra]\n",A.Comm().MyPID());
  ML_Comm* comm;
  ML_Comm_Create(&comm);
  ML_Operator *C_;
  int rv=MatrixMatrix_Multiply(A,comm,&C_);
  Epetra_CrsMatrix_Wrap_ML_Operator(C_,*Comm_,*DomainMap_,C,Copy);
  ML_Operator_Destroy(&C_);
  ML_Comm_Destroy(&comm);
  return rv;
}/*end MatrixMatrix_Multiply*/
#endif
