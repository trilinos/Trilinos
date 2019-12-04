#ifndef HAVE_CONFIG_H
#define HAVE_CONFIG_H
#endif
#include <malloc.h>
#include <cstdio>
#include "ml_config.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_GALERI) && defined(HAVE_ML_AZTECOO)

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Time.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Teuchos_ParameterList.hpp"
#include "ml_MultiLevelPreconditioner.h"
#include "AztecOO.h"

#include "Galeri_Maps.h"
#include "Galeri_Utils.h"
#include "Galeri_CrsMatrices.h"

using namespace Teuchos;
using namespace Galeri;

void PrintLine()
{
  cout << endl;
  for( int i=0 ; i<80 ; ++i )
    cout << "=";
  cout << endl;
  cout << endl;

  return;
}

// ------------------------------------------------------------------------------------------------
// Memory use tracking (gets memory use in MB)
int get_memory() {
  /*** Mallinfo ***/
  struct mallinfo mi = mallinfo();

  printf("Total non-mmapped bytes (arena):       %d\n", mi.arena);
  printf("# of free chunks (ordblks):            %d\n", mi.ordblks);
  printf("# of free fastbin blocks (smblks):     %d\n", mi.smblks);
  printf("# of mapped regions (hblks):           %d\n", mi.hblks);
  printf("Bytes in mapped regions (hblkhd):      %d\n", mi.hblkhd);
  printf("Max. total allocated space (usmblks):  %d\n", mi.usmblks);
  printf("Free bytes held in fastbins (fsmblks): %d\n", mi.fsmblks);
  printf("Total allocated space (uordblks):      %d\n", mi.uordblks);
  printf("Total free space (fordblks):           %d\n", mi.fordblks);
  printf("Topmost releasable block (keepcost):   %d\n", mi.keepcost);
  

  /*** /proc/self/statm ***/
  unsigned long m_vmsize = 0, m_vmrss = 0;
  FILE * fp = fopen( "/proc/self/statm", "r" );
  if (fp)
  { 
    char cbuf[40];
    // run "man proc" to get info on the contents of /proc/[pid]/statm
    fscanf( fp, "%lu %lu %s %s %s %s %s",
            &m_vmsize, &m_vmrss, cbuf, cbuf, cbuf, cbuf, cbuf);

    fclose(fp);
  }
  printf("Virtual memory size (vmsize)       :   %lu\n",m_vmsize);
  printf("Resident set size (vmrss)          :   %lu\n",m_vmrss);

  /*** /prof/self/status ***/
  char vmsize[128];
  char vmpeak[128];
  char vmrss [128];
  char vmhwm[128];
  char line[128];

  FILE *f = fopen("/proc/self/status", "r");
  if (!f) return 0;
  int ct=0;
  
  /* Read memory size data from /proc/pid/status */
  while (ct < 4) {
    if(fgets(line,128,f) == 0)
      break;
    if (!strncmp(line, "VmPeak:", 7)) {
      sprintf(vmpeak,"%s",line); ct++;
    }
    else if (!strncmp(line, "VmSize:", 7)){
      sprintf(vmsize,"%s",line); ct++;
    }
    else if (!strncmp(line, "VmRSS:", 6)) {
      sprintf(vmrss,"%s",line); ct++;
    }
    else if (!strncmp(line, "VmHWM:", 6)) {
      sprintf(vmhwm,"%s",line); ct++;
    }
  }  
  printf("%s", vmpeak);
  printf("%s", vmsize);
  printf("%s", vmrss);
  printf("%s", vmhwm);

  return (mi.hblkhd + mi.usmblks + mi.uordblks)  / 1024;
}


int TestMultiLevelPreconditioner(const char ProblemType[],
				 Teuchos::ParameterList & MLList,
				 Epetra_LinearProblem & Problem, double & TotalErrorResidual,
				 double & TotalErrorExactSol,bool cg=false)
{

  Epetra_MultiVector* lhs = Problem.GetLHS();
  Epetra_MultiVector* rhs = Problem.GetRHS();
  Epetra_RowMatrix* A = Problem.GetMatrix();

  // ======================================== //
  // create a rhs corresponding to lhs or 1's //
  // ======================================== //

  lhs->PutScalar(1.0);
  A->Multiply(false,*lhs,*rhs);

  lhs->PutScalar(0.0);

  Epetra_Time Time(A->Comm());

  // =================== //
  // call ML and AztecOO //
  // =================== //

  AztecOO solver(Problem);

  MLList.set("ML output", 10);

  ML_Epetra::MultiLevelPreconditioner * MLPrec = new ML_Epetra::MultiLevelPreconditioner(*A, MLList, true);

  // tell AztecOO to use this preconditioner, then solve
  solver.SetPrecOperator(MLPrec);

  if(cg) solver.SetAztecOption(AZ_solver, AZ_cg);
  else solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_output, 32);
  solver.SetAztecOption(AZ_kspace, 160);

  solver.Iterate(1550, 1e-12);

  delete MLPrec;

  // ==================================================== //
  // compute difference between exact solution and ML one //
  // ==================================================== //

  double d = 0.0, d_tot = 0.0;

  for( int i=0 ; i<lhs->Map().NumMyElements() ; ++i )
    d += ((*lhs)[0][i] - 1.0) * ((*lhs)[0][i] - 1.0);

  A->Comm().SumAll(&d,&d_tot,1);

  // ================== //
  // compute ||Ax - b|| //
  // ================== //

  double Norm;
  Epetra_Vector Ax(rhs->Map());
  A->Multiply(false, *lhs, Ax);
  Ax.Update(1.0, *rhs, -1.0);
  Ax.Norm2(&Norm);

  string msg = ProblemType;

  if (A->Comm().MyPID() == 0) {
    cout << msg << "......Using " << A->Comm().NumProc() << " processes" << endl;
    cout << msg << "......||A x - b||_2 = " << Norm << endl;
    cout << msg << "......||x_exact - x||_2 = " << sqrt(d_tot) << endl;
    cout << msg << "......Total Time = " << Time.ElapsedTime() << endl;
  }

  TotalErrorExactSol += sqrt(d_tot);
  TotalErrorResidual += Norm;

  return( solver.NumIters() );

}

using namespace Galeri;

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // initialize the random number generator

  int ml_one = 1;
  ML_srandom1(&ml_one);
  // ===================== //
  // create linear problem //
  // ===================== //

  ParameterList GaleriList;
  GaleriList.set("nx", 10);
  GaleriList.set("ny", 10);
  GaleriList.set("nz", 10 * Comm.NumProc());
  GaleriList.set("mx", 1);
  GaleriList.set("my", 1);
  GaleriList.set("mz", Comm.NumProc());

  Epetra_Map* Map = CreateMap("Cartesian3D", Comm, GaleriList);
  Epetra_CrsMatrix* Matrix = CreateCrsMatrix("Laplace3D", Map, GaleriList);
  Epetra_MultiVector* Coords = CreateCartesianCoordinates("3D",Map,GaleriList);
  Epetra_Vector material(Matrix->RowMap());

  Epetra_Vector LHS(*Map);
  Epetra_Vector RHS(*Map);

  Epetra_LinearProblem Problem(Matrix, &LHS, &RHS);

  Teuchos::ParameterList MLList;
  double TotalErrorResidual = 0.0, TotalErrorExactSol = 0.0;
  char mystring[80];

  // ====================== //
  // default options for SA //
  // ====================== //
  if (Comm.MyPID() == 0) PrintLine();

  ML_Epetra::SetDefaults("SA",MLList);
  MLList.set("smoother: type", "Gauss-Seidel");
  strcpy(mystring,"SA");
  TestMultiLevelPreconditioner(mystring, MLList, Problem,
                               TotalErrorResidual, TotalErrorExactSol);


  // ============================================== //
  // default options for SA, efficient symmetric GS //
  // ============================================== //

  if (Comm.MyPID() == 0) PrintLine();

  ML_Epetra::SetDefaults("SA",MLList);
  MLList.set("smoother: type", "Gauss-Seidel");
  MLList.set("smoother: Gauss-Seidel efficient symmetric",true);

  TestMultiLevelPreconditioner(mystring, MLList, Problem,
                               TotalErrorResidual, TotalErrorExactSol,true);

  // ============================== //
  // default options for SA, Jacobi //
  // ============================== //

  if (Comm.MyPID() == 0) PrintLine();

  ML_Epetra::SetDefaults("SA",MLList);
  MLList.set("smoother: type", "Jacobi");

  TestMultiLevelPreconditioner(mystring, MLList, Problem, TotalErrorResidual,
                               TotalErrorExactSol,true);

  // =========================== //
  // default options for SA, Cheby //
  // =========================== //

  if (Comm.MyPID() == 0) PrintLine();

  ML_Epetra::SetDefaults("SA",MLList);
  MLList.set("smoother: type", "Chebyshev");

  TestMultiLevelPreconditioner(mystring, MLList, Problem,
                               TotalErrorResidual, TotalErrorExactSol);



  // =========================== //
  // Specifying Ifpack coarse lists correctly
  // =========================== //
#ifdef HAVE_ML_IFPACK
  if (Comm.MyPID() == 0) PrintLine();
  ML_Epetra::SetDefaults("SA",MLList);

  if(!Comm.MyPID()) {
    MLList.set("ML print initial list",1);
    MLList.set("ML print final list",1);
  }

  MLList.set("smoother: type","ILU");
  MLList.set("coarse: type","ILUT");
  ParameterList &fList = MLList.sublist("smoother: ifpack list");
  fList.set("fact: level-of-fill",1);
  ParameterList &cList = MLList.sublist("coarse: ifpack list");
  cList.set("fact: ilut level-of-fill",1e-2);
  TestMultiLevelPreconditioner(mystring, MLList, Problem,
                               TotalErrorResidual, TotalErrorExactSol);
#endif


  // =========================== //
  // Specifying level sublists
  // =========================== //
  if (Comm.MyPID() == 0) PrintLine();
  ParameterList LevelList;
  ML_Epetra::SetDefaults("SA",LevelList);
  ParameterList &smList = LevelList.sublist("smoother: list (level 0)");
  smList.set("smoother: type","Jacobi");
  smList.set("smoother: sweeps",5);
  ParameterList &smList2 = LevelList.sublist("smoother: list (level 1)");
  smList2.set("smoother: type","symmetric Gauss-Seidel");
  smList2.set("smoother: sweeps",3);
  ParameterList &coarseList = LevelList.sublist("coarse: list");
  coarseList.set("smoother: type","symmetric Gauss-Seidel");
  TestMultiLevelPreconditioner(mystring, LevelList, Problem,
                               TotalErrorResidual, TotalErrorExactSol);

  // =========================== //
  // Ifpack G-S w/ L1
  // =========================== //
#ifdef HAVE_ML_IFPACK
  if (Comm.MyPID() == 0) PrintLine();
  ML_Epetra::SetDefaults("SA",MLList);
  MLList.set("smoother: use l1 Gauss-Seidel",true);
  MLList.set("smoother: type", "Gauss-Seidel");
  TestMultiLevelPreconditioner(mystring, MLList, Problem,
                               TotalErrorResidual, TotalErrorExactSol);
#endif


  // =========================== //
  // Autodetected Line SGS (trivial lines) 
  // =========================== //
  if (Comm.MyPID() == 0) PrintLine();
  ML_Epetra::SetDefaults("SA",MLList);
  MLList.set("smoother: type", "line Gauss-Seidel");
  MLList.set("smoother: line detection threshold",0.1);
  MLList.set("x-coordinates",(*Coords)[0]);
  MLList.set("y-coordinates",(*Coords)[1]);
  MLList.set("z-coordinates",(*Coords)[2]);
  TestMultiLevelPreconditioner(mystring, MLList, Problem,
                               TotalErrorResidual, TotalErrorExactSol);
  

  // =========================== //
  // No QR test                  //
  // =========================== //
  if (Comm.MyPID() == 0) PrintLine();
  ML_Epetra::SetDefaults("SA",MLList);
  MLList.set("aggregation: type","Uncoupled");
  MLList.set("aggregation: do qr",false);
  TestMultiLevelPreconditioner("No QR", MLList, Problem,
                               TotalErrorResidual, TotalErrorExactSol);
  
  // =========================== //
  // Memory Test                 //
  // =========================== //
  int initial_memory = get_memory();
  for (int i=0; i<20; i++) {
    ML_Epetra::SetDefaults("SA",MLList);
    TestMultiLevelPreconditioner(mystring, MLList, Problem,
                                 TotalErrorResidual, TotalErrorExactSol);
  }
  int final_memory = get_memory();
  printf("Memory before = %d after = %d\n",initial_memory,final_memory);

  // =========================== //
  // Material aggregation test
  // =========================== //
  if (Comm.MyPID() == 0) PrintLine();
  material.PutScalar(1.0);
  Teuchos::ParameterList MaterialList;
  ML_Epetra::SetDefaults("SA",MaterialList);
  MaterialList.set("aggregation: material: enable", true);
  MaterialList.set("aggregation: material: threshold", 2.0);
  MaterialList.set("material coordinates", &(material[0]));
  MaterialList.set("x-coordinates",(*Coords)[0]);
  MaterialList.set("y-coordinates",(*Coords)[1]);
  MaterialList.set("z-coordinates",(*Coords)[2]);
  TestMultiLevelPreconditioner("Material", MaterialList, Problem,
                               TotalErrorResidual, TotalErrorExactSol);

  // =========================== //
  // Material & Aux aggregation test
  // =========================== //
  if (Comm.MyPID() == 0) PrintLine();
  material.PutScalar(1.0);
  material[0] = 1000;
  Teuchos::ParameterList MaterialList2;
  ML_Epetra::SetDefaults("SA",MaterialList2);
  MaterialList2.set("aggregation: material: enable", true);
  MaterialList2.set("aggregation: material: threshold", 2.0);
  MaterialList2.set("aggregation: aux: enable", true);
  MaterialList2.set("aggregation: aux: threshold", 1e-10);
  MaterialList2.set("material coordinates", &(material[0]));
  MaterialList2.set("x-coordinates",(*Coords)[0]);
  MaterialList2.set("y-coordinates",(*Coords)[1]);
  MaterialList2.set("z-coordinates",(*Coords)[2]);
  TestMultiLevelPreconditioner("Material", MaterialList2, Problem,
                               TotalErrorResidual, TotalErrorExactSol);


  // =========================== //
  // Classical Test
  // =========================== //
  Teuchos::ParameterList ClassicalList;
  ML_Epetra::SetDefaults("Classical-AMG",ClassicalList);
  TestMultiLevelPreconditioner("Classical-AMG",ClassicalList, Problem,
                               TotalErrorResidual, TotalErrorExactSol);

  // =========================== //
  // Rowsum test (on Heat eqn)   //
  // =========================== //
  if (Comm.MyPID() == 0) PrintLine();
  Epetra_Vector diagonal(Matrix->RowMap());
  Matrix->ExtractDiagonalCopy(diagonal);
  diagonal[0]*=1000;
  Matrix->ReplaceDiagonalValues(diagonal);

  ML_Epetra::SetDefaults("SA",MLList);
  MLList.set("smoother: type", "Chebyshev");
  MLList.set("aggregation: type","Uncoupled");
  MLList.set("aggregation: rowsum threshold", 0.9);

  TestMultiLevelPreconditioner("Rowsum", MLList, Problem,
                               TotalErrorResidual, TotalErrorExactSol);

  // ===================== //
  // print out total error //
  // ===================== //

  if (Comm.MyPID() == 0) {
    cout << endl;
    cout << "......Total error for residual        = " << TotalErrorResidual << endl;
    cout << "......Total error for exact solution  = " << TotalErrorExactSol << endl;
    cout << endl;
  }

  delete Matrix;
  delete Coords;
  delete Map;


  if (TotalErrorResidual > 1e-8) {
    cerr << "Error: `MultiLevelPrecoditioner_Sym.exe' failed!" << endl;
    exit(EXIT_FAILURE);
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  if (Comm.MyPID() == 0)
    cerr << "`MultiLevelPrecoditioner_Sym.exe' passed!" << endl;

  return (EXIT_SUCCESS);
}

#else

#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_MPI
#include "mpi.h"
#endif

int main(int argc, char *argv[])
{
  // still need to deal with MPI, some architecture don't like
  // an exit(0) without MPI_Finalize()
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  puts("Please configure ML with --enable-epetra --enable-teuchos --enable-galeri --enable-aztecoo");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}

#endif /* #if defined(ML_WITH_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_GALERI) */
