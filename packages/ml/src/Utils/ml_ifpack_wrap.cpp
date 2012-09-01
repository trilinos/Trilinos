/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#include "ml_config.h"
#include "ml_include.h"
#if defined(HAVE_ML_IFPACK) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_EPETRA)
#include "ml_utils.h"
#include "ml_epetra.h"
#include "ml_epetra_utils.h"
#include "Epetra_Map.h" 
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h" 
#include "Epetra_VbrMatrix.h" 
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#include "ml_ifpack.h"
#include "ml_ifpack_wrap.h"
#ifdef rst_dump
#include "ml_Ifpack_ML.h"
#endif
// converter from ML_Operator to Epetra_RowMatrix (only wraps)
#include "ml_RowMatrix.h"
// IFPACK factory class
#include "Ifpack.h"
#include "Ifpack_Chebyshev.h"

using namespace ML_Epetra;

static map<void*, bool> MemoryManager;

int ML_Ifpack_Gen(ML *ml, const char* Type, int Overlap, int curr_level, 
                  Teuchos::ParameterList& List, 
                  const Epetra_Comm& Comm, 
                  Ifpack_Handle_Type ** Ifpack_Handle, int reuseSymbolic);

// ====================================================================== 
// MS // This does not work yet with ML_ALL_LEVELS
int ML_Gen_Smoother_Ifpack(ML *ml, const char* Type, int Overlap,
                           int nl, int pre_or_post,
                           void *iList,
                           void *iComm)
{

   int (*fun)(ML_Smoother *, int, double *, int, double *);
   int status = 1;
   char str[80];
   Ifpack_Handle_Type *Ifpack_Handle;
   
   Teuchos::ParameterList List = *((Teuchos::ParameterList *) iList);
   Epetra_Comm *Comm = (Epetra_Comm *) iComm;
#ifdef ML_TIMING
   double         t0;
   t0 = GetClock();
#endif

   fun = ML_Smoother_Ifpack;

   /* Create or reuse IFPACK objects */
   int reuseNumeric  = List.get("reuse numeric factorization",0);
   if (reuseNumeric && ml->pre_smoother[nl].smoother->data != 0) {
     status = 0;
     return status;
   }
   int reuseSymbolic = List.get("reuse symbolic factorization",0);
   if (reuseSymbolic && ml->pre_smoother[nl].smoother->data != 0) {
     Ifpack_Handle = (Ifpack_Handle_Type*) (ml->pre_smoother[nl].smoother->data);
   }
   else {
     Ifpack_Handle = (Ifpack_Handle_Type *) ML_allocate(sizeof(Ifpack_Handle_Type));
     reuseSymbolic = 0;
   } //if (reuseSymbolic)

   status = ML_Ifpack_Gen(ml, Type, Overlap, nl, List, *Comm, &Ifpack_Handle, reuseSymbolic); 
   assert (status == 0); 


   /* This is only used to control the factorization sweeps. I believe */
   /* that when ifpack is used for things like Gauss-Seidel the number */
   /* of sweeps is handled somewhere in IFPACK.                        */

   int sweeps =    List.get("ILU: sweeps", 1);

   /* Sets function pointers */

   if (pre_or_post == ML_PRESMOOTHER) {
     sprintf(str,"IFPACK_pre%d",nl);
     status = ML_Smoother_Set(&(ml->pre_smoother[nl]), (void*)Ifpack_Handle,
			      fun, sweeps, 0.0, str);
     ml->pre_smoother[nl].data_destroy = ML_Smoother_Clean_Ifpack;
#ifdef ML_TIMING
     ml->pre_smoother[nl].build_time = GetClock() - t0;
     ml->timing->total_build_time   += ml->pre_smoother[nl].build_time;
#endif
   }
   else if (pre_or_post == ML_POSTSMOOTHER) {
     sprintf(str,"IFPACK_post%d",nl);
     status = ML_Smoother_Set(&(ml->post_smoother[nl]), 
			      (void*)Ifpack_Handle, fun, sweeps, 0.0, str);
     ml->post_smoother[nl].data_destroy = ML_Smoother_Clean_Ifpack;
#ifdef ML_TIMING
     ml->post_smoother[nl].build_time = GetClock() - t0;
     ml->timing->total_build_time   += ml->post_smoother[nl].build_time;
#endif
   }
   else if (pre_or_post == ML_BOTH) {
     sprintf(str,"IFPACK_pre%d",nl);
     status = ML_Smoother_Set(&(ml->pre_smoother[nl]),
			      (void*)Ifpack_Handle,
			      fun, sweeps,  0.0, str);
     sprintf(str,"IFPACK_post%d",nl);
     status = ML_Smoother_Set(&(ml->post_smoother[nl]),
			      (void*)Ifpack_Handle, fun, sweeps, 0.0, str);
     ml->post_smoother[nl].data_destroy = ML_Smoother_Clean_Ifpack;
#ifdef ML_TIMING
     ml->pre_smoother[nl].build_time = GetClock() - t0;
     ml->timing->total_build_time   += ml->pre_smoother[nl].build_time;
#endif
   }
   else 
     pr_error("ML_Gen_Smoother_Ifpack: unknown pre_or_post choice\n");

   return(status);

}
// ================================================ ====== ==== ==== == =

int ML_Ifpack_Gen(ML *ml, const char* Type, int Overlap, int curr_level, 
                  Teuchos::ParameterList& List, 
                  const Epetra_Comm& Comm, 
                  Ifpack_Handle_Type ** Ifpack_Handle,
                  int reuseSymbolic)
{
# ifdef ML_MPI
  MPI_Comm  ifpackComm;
# else
  int ifpackComm=1;
# endif
  ML_Operator *Ke = &(ml->Amat[curr_level]);
  int hasRows=1;
  bool use_crs=false;

  if (!reuseSymbolic) {
    Epetra_RowMatrix* Ifpack_Matrix;
    Epetra_CrsMatrix* Ifpack_CrsMatrix;
    Epetra_Map*       Ifpack_RowMap;

#ifdef IFPACK_NODE_AWARE_CODE
    // NODE_AWARE always uses use_crs
    use_crs=1;
#endif

    (*Ifpack_Handle)->freeMpiComm = 0;

    if (Ke->type == ML_TYPE_ROW_MATRIX)
    {
      Ifpack_Matrix = (Epetra_RowMatrix*) Ke->data;
      // I have to remember not to delete this guy
      MemoryManager[(void*)Ifpack_Matrix] = false;
    }
    else if(Ke->type == ML_TYPE_CRS_MATRIX)
    {
      Ifpack_Matrix = (Epetra_CrsMatrix*) Ke->data;
      // I have to remember not to delete this guy
      MemoryManager[(void*)Ifpack_Matrix] = false;
    }
    else if(Ke->type == ML_TYPE_VBR_MATRIX)
    {
      Ifpack_Matrix = (Epetra_VbrMatrix*) Ke->data;
      // I have to remember not to delete this guy
      MemoryManager[(void*)Ifpack_Matrix] = false;

    } else {

      // creates the wrapper from ML_Operator to Epetra_RowMatrix
      // (ML_Epetra::RowMatrix). This is a cheap conversion
#     ifdef ML_MPI
      hasRows = MPI_UNDEFINED;
      if (Ke->invec_leng > 0 || Ke->outvec_leng > 0) hasRows = 1;
      MPI_Comm_split(Ke->comm->USR_comm,hasRows,Ke->comm->ML_mypid,&ifpackComm);
      if (hasRows == 1) (*Ifpack_Handle)->freeMpiComm = 1;
#     endif //ifdef ML_MPI
      
      // Do I have CRS storage?
      // Note: MSR keeps everything in cols, so if rowptr isn't null, we're in CRS
      struct ML_CSR_MSRdata* M_= (struct ML_CSR_MSRdata*)ML_Get_MyGetrowData(Ke);
      if(M_ && M_->rowptr) use_crs=true;


      if (hasRows == 1) {
        if(!use_crs){
	// RowMatrix wrapper
          Ifpack_Matrix = new RowMatrix(Ke, 0, false, ifpackComm );
          assert (Ifpack_Matrix != 0);
        }
        else{
          // Uses a CrsMatrix wrapper to enable efficient use of Ifpack on coarser levels
#         ifdef ML_MPI
          Epetra_MpiComm ifpackEpetraComm(ifpackComm);
#         else
          Epetra_SerialComm ifpackEpetraComm;
#         endif
          Ifpack_RowMap = new Epetra_Map(-1,Ke->outvec_leng,0,ifpackEpetraComm);
          Epetra_CrsMatrix_Wrap_ML_Operator(Ke, Ifpack_RowMap->Comm(), *Ifpack_RowMap, &Ifpack_CrsMatrix,View,0);
          assert (Ifpack_CrsMatrix != 0);
          Ifpack_Matrix = Ifpack_CrsMatrix;
          delete Ifpack_RowMap;
        }
        // this guy has to be deleted
        MemoryManager[(void*)Ifpack_Matrix] = true;
      } //if (hasRows == 1)
    } //if (Ke->type == ML_TYPE_ROW_MATRIX) ... else ...

    // we enter the IFPACK world through the factory only
    if (hasRows == 1) {
      Ifpack Factory;
      Ifpack_Preconditioner* Prec;

      // create the preconditioner
      Prec = Factory.Create(Type, Ifpack_Matrix, Overlap);
      Prec->SetParameters(List);
      ML_CHK_ERR(Prec->Compute());
      
      (*Ifpack_Handle)->A_Base = (void *)Prec;
      
      // Grab the lambda's if needed
      if(!strcmp(Type,"Chebyshev")){
        Ifpack_Chebyshev* C=dynamic_cast<Ifpack_Chebyshev*>(Prec);
        assert(C);
        ml->Amat[curr_level].lambda_min=C->GetLambdaMin();
        ml->Amat[curr_level].lambda_max=C->GetLambdaMax();
      }    
    } //if (hasRows==1)
    else {
      (*Ifpack_Handle)->A_Base = 0;
    }

  } else {

    //(*Ifpack_Handle)->A_Base = (void *)Prec;
    Ifpack_Preconditioner* Prec = (Ifpack_Preconditioner*) ((*Ifpack_Handle)->A_Base);
    ML_CHK_ERR(Prec->Compute());

  } //if (!reuseSymbolic) ... else

  
  return 0;
  
} /* ML_Ifpack_Gen */

#ifdef ML_DUMP_IFPACK_FACTORS
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_VectorOut.h"
#include "Ifpack_IC.h"
#include "Ifpack_ICT.h"
#include "Ifpack_ILU.h"
#include "Ifpack_ILUT.h"
#include "Ifpack_AdditiveSchwarz.h"
#endif //ifdef ML_DUMP_IFPACK_FACTORS

// ================================================ ====== ==== ==== == =

int ML_Ifpack_Solve(void * data, double * x, double * rhs )
{
#ifdef rst_dump
   static int stupid_count = 0;
#endif

  Ifpack_Handle_Type *Ifpack_Handle = (Ifpack_Handle_Type *) data;
  if (Ifpack_Handle->A_Base == 0) return 0;

  Ifpack_Preconditioner* Prec = (Ifpack_Preconditioner *)Ifpack_Handle->A_Base;


#ifdef ML_DUMP_IFPACK_FACTORS
  Ifpack_AdditiveSchwarz<Ifpack_IC> *asic = dynamic_cast<Ifpack_AdditiveSchwarz<Ifpack_IC>*>(Prec);
  Ifpack_AdditiveSchwarz<Ifpack_ILU> *asilu = dynamic_cast<Ifpack_AdditiveSchwarz<Ifpack_ILU>*>(Prec);
  Ifpack_AdditiveSchwarz<Ifpack_ICT> *asict = dynamic_cast<Ifpack_AdditiveSchwarz<Ifpack_ICT>*>(Prec);
  Ifpack_AdditiveSchwarz<Ifpack_ILUT> *asilut = dynamic_cast<Ifpack_AdditiveSchwarz<Ifpack_ILUT>*>(Prec);
/*
    const Ifpack_ICT *ict = asict->Inverse();
    ict->Print(cout);
*/
  if (asic != 0) {
    const Ifpack_IC *ic = asic->Inverse();
    const Epetra_Vector &D = ic->D();
    const Epetra_CrsMatrix &U = ic->U();
    EpetraExt::RowMatrixToMatlabFile("IC_Ufactor",U);
    EpetraExt::VectorToMatlabFile("IC_Dfactor",D);
  } else if (asilu != 0) {
    const Ifpack_ILU *ilu = asilu->Inverse();
    const Epetra_CrsMatrix &L = ilu->L();
    const Epetra_Vector &D = ilu->D();
    const Epetra_CrsMatrix &U = ilu->U();
#ifdef rst_dump
    stupid_count++;
    if ((L.RowMatrixRowMap().Comm().MyPID() == 0) && (stupid_count < 6)) {
        Ifpack_ML* MLAndIfpack = (Ifpack_ML *) (Ifpack_Handle->A_Base);
        printf("ILU factors: Nnz(ILU) = %d, Nnz(ILU)/Nnz(A) = %f\n",
               L.NumGlobalNonzeros()+U.NumGlobalNonzeros()-L.NumGlobalRows(),
               ((double)(L.NumGlobalNonzeros()+U.NumGlobalNonzeros()-L.NumGlobalRows()))/
               ((double) MLAndIfpack->Matrix().NumGlobalNonzeros()));
    }
#else
    EpetraExt::RowMatrixToMatlabFile("ILU_Lfactor",L);
    EpetraExt::VectorToMatlabFile("ILU_Dfactor",D);
    EpetraExt::RowMatrixToMatlabFile("ILU_Ufactor",U);
#endif
  } else if (asict != 0) {
    const Ifpack_ICT *ict = asict->Inverse();
    const Epetra_CrsMatrix &H = ict->H();
    EpetraExt::RowMatrixToMatlabFile("ICTfactor",H);
  } else if (asilut != 0) {
    const Ifpack_ILUT *ilut = asilut->Inverse();
    const Epetra_CrsMatrix &L = ilut->L();
    const Epetra_CrsMatrix &U = ilut->U();
    EpetraExt::RowMatrixToMatlabFile("ILUT_Lfactor",L);
    EpetraExt::RowMatrixToMatlabFile("ILUT_Ufactor",U);
  } else {
    if (Prec->Comm().MyPID() == 0) printf("dynamic cast failed!\n");
  }
#ifndef rst_dump
# ifdef HAVE_MPI
  MPI_Finalize();
# endif
  exit(1);
#endif
#endif //ifdef ML_DUMP_IFPACK_FACTORS

  Epetra_Vector Erhs(View, Prec->OperatorRangeMap(), rhs);
  Epetra_Vector Ex(View, Prec->OperatorDomainMap(), x);
  Prec->ApplyInverse(Erhs,Ex); 

  return 0;

} /* ML_Ifpack_Solve */

// ================================================ ====== ==== ==== == =

void ML_Ifpack_Destroy(void * data)
{

  Ifpack_Handle_Type *Ifpack_Handle = (Ifpack_Handle_Type*) data;
  if (Ifpack_Handle->A_Base == 0) {
    ML_free(Ifpack_Handle);
    return;
  }
  //Ifpack_Preconditioner* Prec = (Ifpack_Preconditioner *)Ifpack_Handle;
  Ifpack_Preconditioner* Prec = (Ifpack_Preconditioner*) (Ifpack_Handle->A_Base);

  // a bit nasty, but I don't like the extensive output any more...
  if (ML_Get_PrintLevel() > 10)
    cout << *Prec;

# ifdef ML_MPI
  const Epetra_MpiComm *comm = dynamic_cast<const Epetra_MpiComm*>(&(Prec->Matrix().Comm()));
  if (comm == 0) {
    printf("ML_Ifpack_Destroy: error getting MPI_Comm object\n");
    exit(EXIT_FAILURE);
  }
  MPI_Comm subcomm = comm->GetMpiComm();
# endif

  if (MemoryManager[(void*)(&(Prec->Matrix()))])
  {
    delete &(Prec->Matrix());
    MemoryManager[(void*)(&(Prec->Matrix()))] = false;
  }
  delete Prec;

# ifdef ML_MPI
  if (Ifpack_Handle->freeMpiComm == 1) MPI_Comm_free(&subcomm);
# endif
  ML_free(Ifpack_Handle);

} /* ML_Ifpack_Destroy */

#endif
