#define HAVE_CONFIG_H
#include "Epetra_config.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Time.h"
#include "AztecOO.h"
#include "Epetra_CrsMatrix.h"

#include "Trilinos_Util_ShellOptions.h"
#include "Trilinos_Util_MatrixGallery.h"


int main(int argc, char *argv[])
{
    
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  Epetra_Time Time(Comm);

  Trilinos_Util_CommandLineParser CLP(argc,argv);
  Trilinos_Util_CrsMatrixGallery Gallery("", Comm);

  if( CLP.Has("-problem_type") == false ) CLP.Add("-problem_type", "laplace_2d" ); 
  if( CLP.Has("-problem_size") == false ) CLP.Add("-problem_size", "100" ); 

  Gallery.Set(CLP);

  Epetra_CrsMatrix * Matrix = Gallery.GetMatrix();
  Epetra_LinearProblem * Problem= Gallery.GetLinearProblem();

  AztecOO solver(*Problem);
    
  solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
  solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
  solver.SetAztecOption(AZ_overlap,0);
  solver.SetAztecOption(AZ_subdomain_solve, AZ_lu);
  solver.SetAztecOption(AZ_kspace, 5);
  solver.Iterate(1550, 1e-12);

  double status[AZ_STATUS_SIZE];
  solver.GetAllAztecStatus(status);
  
  // compute the real residual

  double residual, diff;

  Gallery.ComputeResidual(residual);
  Gallery.ComputeDiffBetweenStartingAndExactSolutions(diff);
  
  if( Comm.MyPID()==0 ) 
    cout << "||b-Ax||_2 = " << residual << endl;

  if( Comm.MyPID()==0 ) 
    cout << "||x_exact - x||_2 = " << diff << endl;

  if( Comm.MyPID() == 0 ) 
    cout << "Total Time = " << Time.ElapsedTime() << endl;
#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return 0 ;

}
