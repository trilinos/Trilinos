#include "ml_include.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_GALERI) && defined(HAVE_ML_AZTECOO)

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"
#include "ml_MultiLevelPreconditioner.h"
#include "Epetra_Map.h"
#include "Epetra_LinearProblem.h"

#include <Epetra_Import.h>
#include <Epetra_Export.h>
#include <EpetraExt_CrsMatrixIn.h>
#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_MultiVectorIn.h>
#include <EpetraExt_MultiVectorOut.h>
#include "ml_epetra_utils.h"
#include "ml_RowMatrix.h"


using namespace Teuchos;

EpetraExt::CrsMatrix_SolverMap RowMatrixColMapTrans_;

extern ML_Operator* Epetra_BlkMatRead(char prefix,int, int, ML_Comm *ml_comm,
                                      Epetra_MpiComm &Comm, int , int);
extern  int ML_Smoother_Simple(ML_Smoother *sm,int ,double x[],int , double rhs[]);
extern  int ML_Gen_Smoother_Simple(ML *mlptr, int level);
extern void ML_Smoother_Destroy_Simple(void *data);

struct ML_Simple {
   int nv, np;
   ML *ml00, *ml11;
   ML_Operator *BBt;
   double  *iD;
};
/**************************************************************/
/* A few almost certain changes:                              */
/*     1) ML constructor takes Thrya matrices instead of ML   */
/*        Operators. We might construct an ML operator behind */
/*        the scenes .. so things like ML_Operator and        */
/*        ML_Comm do not appear in the main example.          */
/*     2) ProjectMe() takes and returns Epetra_Row matrices   */
/*     3) ML * arguments to things like ProjectMe() will be   */
/*        changed to ML_Epetra::MultiLevelPreconditioner()    */
/*     4) Smoothers will be invoked by parameter lists.       */
/**************************************************************/

int main(int argc, char *argv[])
{
  int MyPID = 0;
  Epetra_RowMatrix*  A[2];
  int *idummy1, *idummy2; 
  double *ddummy;

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  MyPID = Comm.MyPID();
#else
  Epetra_SerialComm Comm;
#endif
  ML_Comm *ml_comm;
  ML_Comm_Create(&ml_comm);

  /***********************************************/
  /* Read in BlkMatrix, set A[0] to BlkMat(0,0), */
  /* and then read BBt into A[1].                */
  /***********************************************/

  ML_Operator* BlkMat = Epetra_BlkMatRead('s',2, 2, ml_comm, Comm,0,0);
  ML_Operator *tmp = ML_Operator_BlkMatExtract(BlkMat,0,0);

  // CRS or ROW matrix?
  if (Epetra_ML_GetCrsDataptrs(tmp, &ddummy, &idummy1,&idummy2)) 
       A[0] = (Epetra_RowMatrix *) tmp->data;
  else A[0] = dynamic_cast<Epetra_RowMatrix*>((Epetra_CrsMatrix *) tmp->data);

  Epetra_CrsMatrix *ECrsBBt;
  EpetraExt::MatrixMarketFileToCrsMatrix("BBt.mm", Comm, ECrsBBt);
  A[1]= dynamic_cast<Epetra_RowMatrix*>(ECrsBBt);
  A[1]=ML_Epetra::ModifyEpetraMatrixColMap(*(A[1]),RowMatrixColMapTrans_,"BBt");

  /**************************************************/
  /* Build ML parameter lists for subblocks and for */
  /* block 2x2 system.                              */
  /*************************************************/

  ParameterList MLList[2], MLMainList;

  ML_Epetra::SetDefaults("SA",MLMainList);
  ML_Epetra::SetDefaults("SA",MLList[0]);
  ML_Epetra::SetDefaults("SA",MLList[1]);

  MLMainList.set("ML output", 10); 
  MLMainList.set("max levels",3);
  MLMainList.set("smoother: type","do-nothing"); 
  MLMainList.set("coarse: type","do-nothing"); 

  MLList[0].set("ML output", 10);
  MLList[0].set("max levels",3);
  MLList[0].set("PDE equations",2);
  MLList[0].set("smoother: type","do-nothing"); 
  MLList[0].set("coarse: type","do-nothing");

  MLList[1].set("ML output", 10);
  MLList[1].set("max levels",3);
  MLList[1].set("coarse: type","do-nothing"); 
  MLList[1].set("smoother: type","do-nothing"); 

  /*********************************************************/
  /* Constructor for composite AMG. Does AMG on A[k] using */ 
  /* corresponding parameter lists. Then, builds a new AMG */
  /* hierarchy with block diagonal grid transfers. The     */
  /* (k,k)th diagonal block is obtained by extracting the  */
  /* grid transfers associated with AMG on A[k].           */
  /*********************************************************/
  ML_Epetra::MultiLevelPreconditioner* MLPre = 
      new ML_Epetra::MultiLevelPreconditioner(BlkMat, MLMainList, A, MLList, 2);

  /* old-style way of setting smoothers for now */
  void *t=(void *) MLPre->GetML(); ML *mlptr= (ML *) t;
  ML_Gen_Smoother_Simple(mlptr, 0);
  ML_Gen_Smoother_Simple(mlptr, 1);
  ML_Gen_Smoother_Simple(mlptr, 2);

  /* convert BlkMat to an Epetra BlkMat */
  ML_Epetra::RowMatrix EBlkMat(BlkMat,&Comm);
  
  /* read in rhs */
  Epetra_MultiVector *tRHS;
  EpetraExt::MatrixMarketFileToMultiVector("rhs.mm", EBlkMat.Map(), tRHS);
  Epetra_Vector* RHS = dynamic_cast<Epetra_Vector*>(tRHS);

  /* set initial guess */
  Epetra_Vector LHS(EBlkMat.Map()); LHS.PutScalar(0.0);
  for (int i = 0; i < EBlkMat.Map().NumMyElements(); i++) LHS[i] = (*RHS)[i];

  Epetra_LinearProblem Problem(&EBlkMat, &LHS, RHS);
  AztecOO solver(Problem);
  
  solver.SetPrecOperator(MLPre);
  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_output, 1);
  solver.Iterate(24, 1e-10);

  delete MLPre;
  delete A[1];
  ML_Operator_Destroy(&BlkMat);
  ML_Comm_Destroy(&ml_comm);
  ML_spit_it_out();

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}

#else

#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_MPI
#include "mpi.h"
#endif

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  puts("Please configure ML with:");
  puts("--enable-epetra");
  puts("--enable-teuchos");
  puts("--enable-aztecoo");
  puts("--enable-galeri");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  return(EXIT_SUCCESS);
}

#endif
/***********************************************************/
/* Little temporary utility to read in a bunch of Epetra   */
/* CRS matrices and store the whole thing within a block   */
/* ML matrix.                                              */
/***********************************************************/
ML_Operator *Epetra_BlkMatRead(char prefix,int NBlkRows, int NBlkCols, 
          ML_Comm *ml_comm, Epetra_MpiComm &Comm, int MyPID, int DiagOnly) 
{
  ML_Operator* BlkMat = ML_Operator_Create(ml_comm);
  ML_Operator_BlkMatInit(BlkMat, ml_comm, NBlkRows, NBlkCols, 
                         ML_DESTROY_EVERYTHING);
  ML_Operator *tmp;
  Epetra_RowMatrix *EMat;
  Epetra_CrsMatrix *EpetCrs;

  int finfo;
  char str[80];

  for (int iii = 0; iii < NBlkRows ; iii++) {
    for (int jjj = 0; jjj < NBlkCols ; jjj++) {
      if ( (DiagOnly == 0) || (iii == jjj)) {
        sprintf(str,"%c%d_%d.mm",prefix,iii+1,jjj+1);
        finfo = EpetraExt::MatrixMarketFileToCrsMatrix(str, Comm, EpetCrs);
        if (finfo==0) {
          EMat= dynamic_cast<Epetra_RowMatrix*>(EpetCrs);
          EMat = ML_Epetra::ModifyEpetraMatrixColMap(*EMat,RowMatrixColMapTrans_, str);
          tmp = ML_Operator_Create(ml_comm);
          ML_Operator_WrapEpetraMatrix(EMat, tmp);
          ML_Operator_BlkMatInsert(BlkMat, tmp , iii, jjj);
        }
       }
    }
  }
  ML_Operator_BlkMatFinalize(BlkMat);
  return(BlkMat);
}

/********************************************************/
/* Smoother for ML. Must have this proto-type signature */
/* This smoother is motivated by SIMPLE. In particular, */
/* we do the following                                  */
/*       1) x(0) <- smooth( BlkMat(0,0), x(0), rhs(0) ) */
/*       2) r(1) <- BlkMat(1,0)*x(0) - rhs(1)           */
/*       3) x(1) <- smooth( BBt, x(1), r(1) )           */
/*       4) x(0) <- x(0) - iD*BlkMat(0,1)*x(1)          */
/* where                                                */
/*       iD = inv(diag(BlkMat(0,0)))                    */
/********************************************************/
int ML_Smoother_Simple(ML_Smoother *sm,int inlen,double x[],int outlen,
                              double rhs[])
{
   ML_Smoother      *smooth_ptr;
   struct ML_Simple *dataptr;
   ML_Operator      *BlkMat;
   int              nv,np,i;
   double           *iD, *temp, *r;

   smooth_ptr = (ML_Smoother *) sm;
   dataptr    = (struct ML_Simple *) smooth_ptr->smoother->data;
   BlkMat     = smooth_ptr->my_level->Amat;

   nv = dataptr->nv;
   np = dataptr->np;
   iD = dataptr->iD;
   temp = (double *) ML_allocate(sizeof(double)*(nv+np+1));
   r    = (double *) ML_allocate(sizeof(double)*(nv+np+1));

   ML_Smoother_Apply(&(dataptr->ml00->post_smoother[0]),
                  nv, x, nv, rhs, smooth_ptr->init_guess);
   ML_Operator_Apply(ML_Operator_BlkMatExtract(BlkMat,1,0),nv,x,np,temp);

   for (i = 0; i < np; i++) r[i] = temp[i] - rhs[i+nv];

   ML_Smoother_Apply(&(dataptr->ml11->post_smoother[0]),
                     np, &(x[nv]), np, r, ML_NONZERO);

   ML_Operator_Apply(ML_Operator_BlkMatExtract(BlkMat,0,1),np,&(x[nv]),nv,temp);
   for (i = 0; i < nv; i++) x[i] -= (temp[i]*iD[i]);
   
   ML_free(r);
   ML_free(temp);
  return 0;
}
/****************************************************/
/* Simple destroy function invoked by ML destructor */
/****************************************************/
void ML_Smoother_Destroy_Simple(void *data)
{
  struct ML_Simple *widget;

  widget     = (struct ML_Simple *) data;

  ML_free(widget->iD);
  ML_Operator_Destroy(&(widget->BBt));
  ML_Destroy(&(widget->ml00));
  ML_Destroy(&(widget->ml11));
  ML_free(widget);

}
/****************************************************/
/* Construct an ML smoother object corresponding to */
/* simple.                                          */
/****************************************************/
int ML_Gen_Smoother_Simple(ML *mlptr, int level)
{
  ML_Operator *Fmat, *DinvBt, *BBt, *BlkMat;
  double *Dinv, *val;
  int    allocated = 100, *bindx, row_length;
  struct ML_Simple *widget;
  ML_Smoother *smooth_ptr;
   
  BlkMat = &(mlptr->Amat[level]);
  Fmat = ML_Operator_BlkMatExtract(BlkMat,0,0);


  /* construct Dinv */

  Dinv  = (double  *) ML_allocate(sizeof(double)*(Fmat->outvec_leng+1));
  val   = (double  *) ML_allocate(sizeof(double)*allocated);
  bindx = (int     *) ML_allocate(sizeof(int   )*allocated);
  for (int row = 0; row < Fmat->outvec_leng; row++) {
    Dinv[row] = 1.;
    ML_get_matrix_row(Fmat, 1, &row, &allocated, &bindx,&val,&row_length,0);
    for (int j = 0; j < row_length; j++)
       if ( (bindx[j] == row) && (val[j] != 0.)) Dinv[row] = 1./val[j];
  }
  ML_free(bindx); ML_free(val);

  /* Build  B*Dinv*B'. On finest level, this requires multiplying */
  /* the three matrices together. On coarse levels this is done   */
  /* by projecting the fine level BBt.                            */

  if (level == 0) {
    BBt = ML_Operator_Create(Fmat->comm);
    DinvBt=ML_Operator_ImplicitlyVScale(ML_Operator_BlkMatExtract(BlkMat,0,1),Dinv,0);
    ML_2matmult(ML_Operator_BlkMatExtract(BlkMat,1,0),DinvBt,BBt,ML_MSR_MATRIX);
    ML_Operator_Destroy(&DinvBt);
  }
  else {
    smooth_ptr= &(mlptr->post_smoother[level-1]);
    widget    = (struct ML_Simple *) smooth_ptr->smoother->data;
    BBt       = ProjectMe(mlptr, 1, level-1, widget->BBt, ML_MSR_MATRIX);
  }

  /* Make smoother widget. ML is a bit funny in that the     */
  /* smoothers for the SIMPLE submatrices must be associated */
  /* with a multigrid hierarchy. Thus, we create a 1-level   */
  /* AMG hierarchy to hold the smoothers for the (0,0) and   */
  /* the (1,1) block.                                        */
  /* NOTE: if this is the coarsest level we use Amesos KLU   */
  /*       instead of Gauss_Seidel as the smoother.          */

  ML *ml00, *ml11;
  ML_Create(&ml00,1);
  ML_Operator_halfClone_Init( &(ml00->Amat[0]), Fmat);
  if (mlptr->Rmat[level].to != NULL) 
    ML_Gen_Smoother_GaussSeidel(ml00, 0, ML_POSTSMOOTHER,1,  1.);
  else ML_Gen_Smoother_Amesos(  ml00, 0,   ML_AMESOS_KLU,1, 0.0);

  ML_Create(&ml11,1);
  ML_Operator_halfClone_Init( &(ml11->Amat[0]), BBt);
  if (mlptr->Rmat[level].to != NULL) 
    ML_Gen_Smoother_GaussSeidel(ml11, 0, ML_POSTSMOOTHER,1, 1.);
  else ML_Gen_Smoother_Amesos(  ml11, 0,   ML_AMESOS_KLU,1, 0.0);

  widget = (struct ML_Simple *) ML_allocate(sizeof(struct ML_Simple));
  widget->nv   = Fmat->invec_leng;
  widget->np   = BBt->invec_leng;
  widget->ml00 = ml00;
  widget->ml11 = ml11;
  widget->BBt  = BBt;
  widget->iD   = Dinv;
    
  if (mlptr->Rmat[level].to != NULL) 
    ML_Smoother_Set(&(mlptr->pre_smoother[level]), (void *) widget, 
                    ML_Smoother_Simple, 1, 1.0, "hi");

  ML_Smoother_Set(&(mlptr->post_smoother[level]), (void *) widget, 
                  ML_Smoother_Simple, 1, 1.0, "hi");
  mlptr->post_smoother[level].data_destroy = ML_Smoother_Destroy_Simple;
}
