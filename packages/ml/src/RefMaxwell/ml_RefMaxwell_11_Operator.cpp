#include "ml_config.h"
#include "ml_epetra.h"
#include "ml_RefMaxwell_11_Operator.h"
#include "ml_epetra_utils.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_EPETRAEXT)
#include "EpetraExt_Transpose_RowMatrix.h"
#include "EpetraExt_SolverMap_CrsMatrix.h"
#include "EpetraExt_MatrixMatrix.h" //haq


#define BASE_IDX 0
#define NO_OUTPUT

#ifdef NO_OUTPUT
#define MVOUT2(x,y,z) ;
#define MVOUT(x,y) ;
#define IVOUT(x,y) ;
#define Epetra_CrsMatrix_Print(x,y) ;
#define ML_Matrix_Print(w,x,y,z) ;
#else
extern void Epetra_CrsMatrix_Print(const Epetra_CrsMatrix& A, char* of);
extern void MVOUT (const Epetra_MultiVector & A, char * of);
extern void IVOUT(const Epetra_IntVector & A, char * of);
extern void MVOUT2(const Epetra_MultiVector & A,char* pref,int idx);
extern void ML_Matrix_Print(ML_Operator *ML,const Epetra_Comm &Comm,const Epetra_Map &Map, char *fname);
#endif


extern void print_stats(const Epetra_CrsMatrix& A, char *label);//haq
// ================================================ ====== ==== ==== == = 
// Constructor
ML_Epetra::ML_RefMaxwell_11_Operator::ML_RefMaxwell_11_Operator(const Epetra_CrsMatrix& SM_Matrix,    //S+M
                                                             const Epetra_CrsMatrix& D0_Matrix,    //T or D0
                                                             const Epetra_CrsMatrix& M0inv_Matrix, //M0^{-1}
                                                             const Epetra_CrsMatrix& M1_Matrix):   //M1(1)
  SM_Matrix_(&SM_Matrix),Addon_Matrix_(0)
#ifdef USE_CORE_MATRIX
  ,Core_Matrix_(0),Core_Matrix_Reindex_(0),Core_Matrix_Reindexer_(0)
#else
  ,D0T_Matrix_(0)
#endif
{
  Label_=strdup("ML_RefMaxwell_11_Operator");
  Comm_ = &(SM_Matrix_->Comm());
  DomainMap_ = &(SM_Matrix_->OperatorDomainMap());
  RangeMap_ = &(SM_Matrix_->OperatorRangeMap());


#ifdef USE_CORE_MATRIX  
  printf("Building Core Matrix\n");
  /* Build the DO M0 D0^T Core */
  //  Epetra_CrsMatrix *core_temp[2] = {(Epetra_CrsMatrix*)&D0_Matrix,(Epetra_CrsMatrix*)&M0inv_Matrix};
  //  Epetra_Multi_CrsMatrix CoreTemp(2,core_temp);
  //  int rv=CoreTemp.MatrixMatrix_Multiply(D0_Matrix_transpose,&Core_Matrix_);a
  //  if(rv) printf("[%d] ERROR: Core Matrix formation failed with error code %d\n",Comm_->MyPID(),rv);
  //  else  printf("[%d] ML_RefMaxwell_11_Operator: Core Built\n",Comm_->MyPID());
  
  /* Build the D0 M0 D0^T Core with EpetraExt */
  Epetra_CrsMatrix temp1(Copy,M0inv_Matrix.RangeMap(),BASE_IDX);
  Core_Matrix_=new Epetra_CrsMatrix(Copy,*RangeMap_,BASE_IDX);
  
  EpetraExt::MatrixMatrix::Multiply(M0inv_Matrix,false,D0_Matrix,true,temp1);
  EpetraExt::MatrixMatrix::Multiply(D0_Matrix,false,temp1,false,*Core_Matrix_);
  Core_Matrix_->OptimizeStorage();
  //NTS: For whatever reason ML's wrap routines die *horribly* when we play with
  //the transposed matrix created by EpetraExt::RowMatrix_Transpose.  This
  //should be fixed eventually.  As for now, we'll just use EpetraExt to build
  //the core.

  // Deal w/ reindexing garbage for the core matrix
  Epetra_Map RangeMap_Consec(-1,DomainMap_->NumMyElements(),BASE_IDX,*Comm_);
  Core_Matrix_Reindexer_=new EpetraExt::CrsMatrix_Reindex(RangeMap_Consec);    
  Core_Matrix_Reindex_=&((*Core_Matrix_Reindexer_)(*Core_Matrix_));

#ifndef NO_OUTPUT
  Epetra_CrsMatrix_Print(*Core_Matrix_,"core_matrix.dat");
#endif
  
  /* Build the Epetra_Multi_CrsMatrix */
  Addon_Matrix_=new Epetra_CrsMatrix*[3];
  Addon_Matrix_[0]=Addon_Matrix_[2]=(Epetra_CrsMatrix*)&M1_Matrix;
  Addon_Matrix_[1]=Core_Matrix_Reindex_;  
  Addon_=new Epetra_Multi_CrsMatrix(3,Addon_Matrix_);
#else

  /* Transpose D0 */
  D0_Matrix_Transposer_= new EpetraExt::RowMatrix_Transpose(true,(Epetra_Map*)&M0inv_Matrix.OperatorRangeMap());
  D0T_Matrix_= dynamic_cast<Epetra_CrsMatrix*>( & ((*D0_Matrix_Transposer_)((Epetra_CrsMatrix&)D0_Matrix)));
  D0T_Matrix_= dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(*D0T_Matrix_,D0T_Matrix_Trans_,"D0T",false));
  
  /* Build the Epetra_Multi_CrsMatrix */
  Addon_Matrix_=new Epetra_CrsMatrix*[5];
  Addon_Matrix_[0]=Addon_Matrix_[4]=(Epetra_CrsMatrix*)&M1_Matrix;
  Addon_Matrix_[1]=(Epetra_CrsMatrix*)&D0_Matrix;Addon_Matrix_[3]=D0T_Matrix_;
  Addon_Matrix_[2]=(Epetra_CrsMatrix*)&M0inv_Matrix;
  Addon_=new Epetra_Multi_CrsMatrix(5,Addon_Matrix_);
#endif  
}/*end Constructor*/ 

// ================================================ ====== ==== ==== == = 
// Destructor
ML_Epetra::ML_RefMaxwell_11_Operator::~ML_RefMaxwell_11_Operator()
{
  if(Label_) free(Label_);
  if(Addon_Matrix_) delete [] Addon_Matrix_;
  if(Addon_) delete Addon_;
#ifdef USE_CORE_MATRIX
  if(Core_Matrix_) delete Core_Matrix_;
  if(Core_Matrix_Reindexer_) delete Core_Matrix_Reindexer_;
#else
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
// Computes C= <me> * A
int  ML_Epetra::ML_RefMaxwell_11_Operator::MatrixMatrix_Multiply(const Epetra_CrsMatrix & A, ML_Comm *comm, ML_Operator **C) const
{  
  ML_Operator *SM_ML,*A_ML,*temp1,*temp2;

  /* General Stuff */  
  ML_Comm* temp = global_comm;
  A_ML  = ML_Operator_Create(comm);
  *C = ML_Operator_Create(comm);
  //  ML_Operator_WrapEpetraCrsMatrix((Epetra_CrsMatrix*)&A,A_ML);  
  ML_Operator_WrapEpetraMatrix((Epetra_CrsMatrix*)&A,A_ML);  

  
  /* Do the SM part */
  SM_ML = ML_Operator_Create(comm);
  temp1 = ML_Operator_Create(comm);
  //  ML_Operator_WrapEpetraCrsMatrix((Epetra_CrsMatrix*)SM_Matrix_,SM_ML);
  ML_Operator_WrapEpetraMatrix((Epetra_CrsMatrix*)SM_Matrix_,SM_ML);
  ML_2matmult(SM_ML,A_ML,temp1,ML_CSR_MATRIX);
  ML_Matrix_Print(temp1,*Comm_,*RangeMap_,"smp.dat");
    
  //#define SANITY_CHECK
#ifdef SANITY_CHECK
  /* DEBUG */
  Epetra_CrsMatrix C_EP3(Copy,*DomainMap_,BASE_IDX);
  EpetraExt::MatrixMatrix::Multiply(*SM_Matrix_,false,A,false,C_EP3);
  Epetra_CrsMatrix diff(C_EP3), *C_EP2;
  Epetra_CrsMatrix_Wrap_ML_Operator(temp1,*Comm_,*DomainMap_,&C_EP2,Copy,SM_Matrix_->IndexBase());
  C_EP2->OptimizeStorage();
  EpetraExt::MatrixMatrix::Add(*C_EP2,false,1.0,diff,-1.0);

  
  double norm1=C_EP3.NormInf();  
  double norm2=C_EP2->NormInf();
  double norm3=diff.NormInf()/norm1;
  if(A.Comm().MyPID()==0)
     printf("RM11: diff = %6.4e (%6.4e/%6.4e)\n",norm3,norm1,norm2);
  delete C_EP2;
  /*end DEBUG*/
#endif

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
  ML_Operator *C_;
  int rv=MatrixMatrix_Multiply(A,comm,&C_);
  Epetra_CrsMatrix_Wrap_ML_Operator(C_,*Comm_,*DomainMap_,C,Copy,A.IndexBase());
  ML_Operator_Destroy(&C_);
  ML_Comm_Destroy(&comm);
  return rv;
}/*end MatrixMatrix_Multiply*/
#endif
