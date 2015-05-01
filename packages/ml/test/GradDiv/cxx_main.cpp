/*
 * Goal of this test:
 * - verify that the GradDiv Solver functions correctly
 *
 * This test will:
 * - Read matrices from disk
 * - Form a GradDiv preconditioner
 * - Test several combinations of options
 *
 * \date Sept. 22, 2011
 *
 * \author Chris Siefert, 01423
 *
 */

#include <iostream>
#include <math.h>
#ifndef ICL
#ifndef WIN32
#include <unistd.h>
#endif
#endif
#include "ml_include.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_EPETRAEXT)


#include "EpetraExt_MatrixMatrix.h"
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_SolverMap_CrsMatrix.h"
#include "EpetraExt_MultiVectorIn.h"
#include "ml_epetra_utils.h"

#include "ml_GradDiv.h"
#include "Epetra_LinearProblem.h"
#include "Teuchos_ParameterList.hpp"
#include "AztecOO.h"

using namespace ML_Epetra;

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#define Epetra_ActiveComm Epetra_MpiComm
#else
#include "Epetra_SerialComm.h"
#define Epetra_ActiveComm Epetra_SerialComm
#endif
/*******************************************************/
int MatlabFileToMultiVector( const char *filename, const Epetra_BlockMap & map, int N, Epetra_MultiVector * & A);
int dim=3;



/*******************************************************/
Teuchos::ParameterList Build_Teuchos_List(double *coordx, double*coordy, double* coordz,const char *str_tag, const char* str_val,const char *int_tag, int int_val){
  Teuchos::ParameterList List_Coarse, RMList;

  /* Pass in given options */
  List_Coarse.set("PDE equations",dim);
  List_Coarse.set("x-coordinates",coordx);
  List_Coarse.set("y-coordinates",coordy);
  List_Coarse.set("z-coordinates",coordz);
  List_Coarse.set("ML output",10);

  Teuchos::ParameterList List11,List11c,List22,List22c;
  ML_Epetra::UpdateList(List_Coarse,List11,true);
  ML_Epetra::UpdateList(List_Coarse,List22,true);
  ML_Epetra::UpdateList(List_Coarse,List11c,true);
  ML_Epetra::UpdateList(List_Coarse,List22c,true);

  if(int_tag) List11c.set(int_tag,int_val);
  if(str_tag) List11c.set(str_tag,str_val);
  if(int_tag) List22c.set(int_tag,int_val);
  if(str_tag) List22c.set(str_tag,str_val);


  /* Setup basic list structure */
  List11.set("face matrix free: coarse",List11c);
  List22.set("edge matrix free: coarse",List22c);
  RMList.setName("graddiv list");
  RMList.set("graddiv: 11list",List11);
  RMList.set("graddiv: 22list",List22);
  if(int_tag) RMList.set(int_tag,int_val);
  if(str_tag) RMList.set(str_tag,str_val);

  return RMList;
}


/*******************************************************/
void gd_test_additive(Epetra_ActiveComm &Comm,
                       Teuchos::ParameterList & List,
                       Epetra_CrsMatrix &K2,
                       const Epetra_CrsMatrix &FaceNode,
                       const Epetra_CrsMatrix &D1,
		       const Epetra_CrsMatrix &D0,
  		       const Epetra_CrsMatrix &K0,
                       const Epetra_Vector &x_exact,
                       const Epetra_Vector &x0,
                       const Epetra_Vector &b){

  cout<<"*************************************"<<endl;
  cout<<List<<endl;
  cout<<"*************************************"<<endl;

  GradDivPreconditioner PrecRF(K2,FaceNode,D1,D0,K0,List);
  /*Teuchos::ParameterList List2;
  SetDefaults("SA",List2);
  List2.set("ML output",10);
  MultiLevelPreconditioner PrecRF(K2,List2);*/

  Epetra_Vector x0_(x0);

  Epetra_LinearProblem Problem(&K2, &x0_, (Epetra_MultiVector*)&b);
  AztecOO solver(Problem);
  solver.SetPrecOperator(&PrecRF);

  solver.SetAztecOption(AZ_solver, AZ_cg);
  solver.SetAztecOption(AZ_conv,AZ_r0);
  solver.SetAztecOption(AZ_output,10);
  solver.Iterate(100,1e-8);

  /* Check out the solution */
  double nxe,nd;
  x_exact.Norm2(&nxe);
  Epetra_Vector diff(x_exact);
  diff.Update(1.0,x0_,-1.0);
  diff.Norm2(&nd);
  if(Comm.MyPID()==0) printf("||sol-exact||/||exact||=%6.4e\n",nd/nxe);
  if(nd/nxe > 1e-6) exit(1);


}



/*******************************************************/
void matrix_read(Epetra_ActiveComm &Comm){
  Epetra_CrsMatrix *K2e,*FaceNodee,*D1e,*D0e,*K0e;

  /* Read Matrices */
  EpetraExt::MatlabFileToCrsMatrix("k2.dat" ,Comm,K2e);
  EpetraExt::MatlabFileToCrsMatrix("facenode.dat",Comm,FaceNodee);
  EpetraExt::MatlabFileToCrsMatrix("d1.dat",Comm,D1e);
  EpetraExt::MatlabFileToCrsMatrix("d0.dat",Comm,D0e);
  EpetraExt::MatlabFileToCrsMatrix("k0.dat" ,Comm,K0e);

  /* Optimize Storage */
  K2e->OptimizeStorage();
  FaceNodee->OptimizeStorage();
  D1e->OptimizeStorage();
  D0e->OptimizeStorage();
  K0e->OptimizeStorage();

  /* Fix Column Maps */
  EpetraExt::CrsMatrix_SolverMap K2_CMT, FaceNode_CMT, D1_CMT, D0_CMT,K0_CMT;
  Epetra_CrsMatrix *K2,*FaceNode,*D1,*D0,*K0;
  K2=dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(*K2e, K2_CMT ));
  FaceNode=dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(*FaceNodee,FaceNode_CMT ));
  D1=dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(*D1e,D1_CMT));
  D0=dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(*D0e,D0_CMT));
  K0=dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(*K0e,K0_CMT));

  /* Optimize Storage */
  K2->OptimizeStorage();
  FaceNode->OptimizeStorage();
  D1->OptimizeStorage();
  D0->OptimizeStorage();
  K0->OptimizeStorage();

  printf("size(K2)=%d size(D1)=%dx%d size(D0)=%dx%d size(K0)=%d\n",K2->NumGlobalRows(),D1->NumGlobalRows(),D0->NumGlobalRows(),
	 D0->NumGlobalCols(),D1->NumGlobalCols(),K0->NumGlobalRows());

  /* Build RHS */
  const Epetra_Map &FaceMap=K2->DomainMap();
  const Epetra_Map &NodeMap=K0->DomainMap();
  Epetra_Vector rhs(FaceMap,false);
  Epetra_Vector x_exact(FaceMap,false);
  x_exact.PutScalar(1.0);
  K2->Apply(x_exact,rhs);

  /* Read in coordinates*/
  Epetra_MultiVector *coordx,*coordy,*coordz;
  MatlabFileToMultiVector("coordx.dat",NodeMap,dim,coordx);
  MatlabFileToMultiVector("coordy.dat",NodeMap,dim,coordy);
  MatlabFileToMultiVector("coordz.dat",NodeMap,dim,coordz);

  /* Build Lists */
  Teuchos::ParameterList List_2level = Build_Teuchos_List(coordx->Values(),coordy->Values(),coordz->Values(),"coarse: type","Amesos-KLU","max levels",1);
  Teuchos::ParameterList List_Cheby  = Build_Teuchos_List(coordx->Values(),coordy->Values(),coordz->Values(),0,0,0,0);

  /* Do Tests */
  Epetra_Vector lhs(FaceMap,true);

  gd_test_additive(Comm,List_2level,*K2,*FaceNode,*D1,*D0,*K0,x_exact,lhs,rhs);

  lhs.PutScalar(0.0);
  gd_test_additive(Comm,List_Cheby,*K2,*FaceNode,*D1,*D0,*K0,x_exact,lhs,rhs);

  delete K2; if(K2e!=K2) delete K2e;
  delete FaceNode; if(FaceNode!=FaceNodee) delete FaceNodee;
  delete D1; if(D1!=D1e) delete D1e;
  delete D0; if(D0!=D0e) delete D0e;
  delete K0; if(K0!=K0e) delete K0e;
  delete coordx; delete coordy; delete coordz;

}


/*******************************************************/
int MatlabFileToMultiVector( const char *filename, const Epetra_BlockMap & map, int N, Epetra_MultiVector * & A) {

  const int lineLength = 1025;
  char line[lineLength];

  FILE * handle = 0;

  handle = fopen(filename,"r");  // Open file
  if (handle == 0)
    EPETRA_CHK_ERR(-1); // file not found

  // Next, strip off header lines (which start with "%")
  //  do {
  //    if(fgets(line, lineLength, handle)==0) return(-1);
  //  } while (line[0] == '%');

  // Compute the offset for each processor for when it should start storing values
  int numMyPoints = map.NumMyPoints();
  int offset=0;
  map.Comm().ScanSum(&numMyPoints, &offset, 1); // ScanSum will compute offsets for us
  offset -= numMyPoints; // readjust for my PE
  if(map.Comm().NumProc() == 1) offset=0;//CMS

  // Now construct vector/multivector
  if (N==1)
    A = new Epetra_Vector(map);
  else
    A = new Epetra_MultiVector(map, N);

  double ** Ap = A->Pointers();

  // Now read in lines that we will discard
  for (int i=0; i<offset; i++)
    if(fgets(line, lineLength, handle)==0) return(-3);
  for (int i=0; i<numMyPoints; i++) {
    for (int j=0; j<N; j++) {
      double * v = Ap[j];
      // Now read in each value and store to the local portion of the the  if the row is owned.
      double V;
      if(fscanf(handle, "%le ", &V)==0) return(-5);
      v[i] = V;
    }
  }

  if (fclose(handle)) return(-1);

  return(0);
}


/*******************************************************/
int main(int argc, char* argv[]){
  /* Initialize */
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  matrix_read(Comm);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}/*end main*/





#else

#ifdef HAVE_MPI
#include "mpi.h"
#endif

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  puts("This test requires:");
  puts("--enable-epetra");
  puts("--enable-aztecoo");
  puts("--enable-epetraext");
  puts("--enable-ifpack");
  puts("--enable-teuchos");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  // not to break tests
  return(EXIT_SUCCESS);
}

#endif
