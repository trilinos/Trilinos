/*!
 *  \file ml_MultiLevelPreconditioner_Analyze.cpp
 *
 *  \brief Visualization utilities for MultiLevelPreconditioner class
 *
 *  \author Marzio Sala, SNL, 9214
 *
 *  \date Last update do Doxygen: 09-Aug-04
 *
 */

#include "ml_common.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Epetra_Time.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_FECrsMatrix.h"
#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "ml_epetra_preconditioner.h"
#ifdef HAVE_ML_AZTECOO
#include "AztecOO.h"
#endif
#include <iomanip>

// ============================================================================
static void MLP_print(int count, char * str, double status[AZ_STATUS_SIZE], double time)
{

  cout << "#" << count;
  if( count < 10  )      cout << ".....";
  else if( count < 100 ) cout << "....";
  else cout << "...";

  cout.width(30); cout.setf(ios::left); cout.fill('.');
  cout << str;
  cout.width(10); cout.setf(ios::left); cout.fill('.');
  cout << (int)status[AZ_its];
  cout.width(15); cout.setf(ios::left); cout.fill('.');
  cout << status[AZ_scaled_r];
  cout.width(15); cout.setf(ios::left); cout.fill('.');
  cout << time;
  
  if( status[AZ_why] == AZ_normal         ) cout << "N";
  else if( status[AZ_why] == AZ_maxits    ) cout << "M";
  else if( status[AZ_why] == AZ_loss      ) cout << "L";
  else if( status[AZ_why] == AZ_ill_cond  ) cout << "I";
  else if( status[AZ_why] == AZ_breakdown ) cout << "B";

  cout << endl;
}

// ============================================================================
int ML_Epetra::MultiLevelPreconditioner::AnalyzeMatrix(char * Defaults, bool IsSymmetric)
{

  if( RowMatrix_ == 0 ) return -1; // Matrix not yet set

  bool Properties  = List_.get("analysis: matrix properties", true);
  bool Eigenvalues = List_.get("analysis: matrix eigenvalues", false);
  bool Smoother    = List_.get("analysis: smoothers", true);
  // bool Coarsening  = List_.get("analysis: coarsening", true);

  // =============================== //
  // few analysis on matrix property //
  // =============================== //

  if( Properties ) AnalyzeMatrixProperties(Defaults,IsSymmetric);

  // ======================== //
  // compute some eigenvalues //
  // ======================== //

  if( Eigenvalues ) AnalyzeMatrixEigenvalues(Defaults,IsSymmetric);
  // 1.- smallest element
  // 2.- largest element
  // 3.- is the matrix positive/negative?
  // 4.- simple estimation of SPD
  // 5.- computation of eigenvalues using Anasazi
  // 6.- scaling

  // ============================= //
  // effect of different smoothers //
  // ============================= //

  if( Smoother ) AnalyzeMatrixSmoothers(Defaults,IsSymmetric);

  // ============================== //
  // effect of different coarsening //
  // ============================== //
  // I am not sure of this right now
  // if( Coarsening ) AnalyzeMatrixCoarsening(IsSymmetric);

  return 0;
}
// ============================================================================
void ML_Epetra::MultiLevelPreconditioner::AnalyzeMatrixProperties(char * Defaults,
								  bool IsSymmetric)
{
  return;
}
// ============================================================================
void ML_Epetra::MultiLevelPreconditioner::AnalyzeMatrixEigenvalues(char * Defaults, 
								   bool IsSymmetric)
{
  return;
}
// ============================================================================
void ML_Epetra::MultiLevelPreconditioner::AnalyzeMatrixSmoothers(char * Defaults,
								 bool IsSymmetric)
{
  
  int MaxIters = List_.get("analysis: max iters",500);
  double Tol   = List_.get("analysis: tolerance", 1e-5);
  double status[AZ_STATUS_SIZE];
  char smoother[80];
  int sweeps = 1;
  Epetra_Time Time(Comm());
  int count=0;
  double ReqTime, BestIters = 1000000.0, BestTime = 1000000.0;
  int BestItersCount = -1, BestTimeCount = -1;
  
  // create a new MultiLevePreconditioner based on the same matrix
  ML_Epetra::MultiLevelPreconditioner * yo;

  // ========================== //
  // create the AztecOO problem //
  // ========================== //
  Epetra_Vector LHS(Map());
  Epetra_Vector RHS(Map());

  Epetra_LinearProblem Problem(const_cast<Epetra_RowMatrix*>(RowMatrix_),&LHS, &RHS);

  AztecOO solver(Problem);

  if( IsSymmetric ) solver.SetAztecOption(AZ_solver, AZ_cg);
  else              solver.SetAztecOption(AZ_solver, AZ_gmres); 
  solver.SetAztecOption(AZ_kspace, 50);
  solver.SetAztecOption(AZ_conv, AZ_r0);
  solver.SetAztecOption(AZ_output, AZ_none);

  // output
 
  if( Comm().MyPID() == 0 ) {
    cout << "@@@" << endl;
    cout << "@@@ Analysis of ML parameters (smoothers)" << endl;
    cout << "@@@" << endl << endl;;
    cout << "@@@ maximum iterations = " << MaxIters << endl;
    cout << "@@@ tolerance          = " << Tol << endl;
    cout << "@@@ Using default options";
    if( Defaults ) cout << " (" << Defaults << ")";
    cout << endl;
    cout << endl;
    cout << "count  ";
    cout.width(30); cout.setf(ios::left); cout.fill('.');
    cout << "smoother type";
    cout.width(10); cout.setf(ios::left); cout.fill('.');
    cout << "its";
    cout.width(15); cout.setf(ios::left); cout.fill('.');
    cout << "||r||/||r_0||";
    cout.width(15); cout.setf(ios::left); cout.fill('.');
    cout << "time (s)" << endl;
  }

  // ====== //
  // Jacobi //
  // ====== //

  if( Comm().MyPID() == 0 ) cout << endl << "- Jacobi" << endl;

  for( double omega=0.25 ; omega<1.5 ; omega+=0.25)
  {
    Time.ResetStartTime();

    Teuchos::ParameterList NewList;
    ML_Epetra::SetDefaults(Defaults, NewList);
    NewList.set("output", 0);
    NewList.set("smoother: type", "Jacobi");
    NewList.set("smoother: damping", omega);
    NewList.set("smoother: sweeps", sweeps);
    yo = new ML_Epetra::MultiLevelPreconditioner(*RowMatrix_,NewList, true);
    assert( yo != 0 );
    solver.SetPrecOperator(yo);

    LHS.PutScalar(0.0);
    RHS.Random();

    solver.Iterate(MaxIters,Tol);
    solver.GetAllAztecStatus(status);
    sprintf(smoother,"n=%d, omega=%5.2e", sweeps, omega);
    ReqTime = Time.ElapsedTime();
    if( ReqTime < BestTime ) BestTimeCount = count;
    if( (int) status[AZ_its] < BestIters ) BestItersCount = count;
    if( Comm().MyPID() == 0 ) MLP_print(count++,smoother,status,ReqTime);

    delete yo;
  }
  
  // ============ //
  // Gauss-Seidel //
  // ============ //

  if( Comm().MyPID() == 0 ) cout << endl << "- Gauss-Seidel" << endl;

  for( double omega=0.25 ; omega<1.5 ; omega+=0.25)
  {

    Time.ResetStartTime();

    Teuchos::ParameterList NewList(List_);
    ML_Epetra::SetDefaults(Defaults, NewList);
    NewList.set("output", 0);

    NewList.set("smoother: type", "Gauss-Seidel");
    NewList.set("smoother: damping", omega);
    NewList.set("smoother: sweeps", sweeps);
    yo = new ML_Epetra::MultiLevelPreconditioner(*RowMatrix_,NewList, true);
    assert( yo != 0 );
    solver.SetPrecOperator(yo);

    LHS.PutScalar(0.0);
    RHS.Random();

    solver.Iterate(MaxIters,Tol);
    solver.GetAllAztecStatus(status);
    sprintf(smoother,"n=%d, omega=%5.2e", sweeps, omega);
    ReqTime = Time.ElapsedTime();
    if( ReqTime < BestTime ) BestTimeCount = count;
    if( (int) status[AZ_its] < BestIters ) BestItersCount = count;
    if( Comm().MyPID() == 0 ) MLP_print(count++,smoother,status,ReqTime);

    delete yo;
  }

  // ====================== //
  // symmetric Gauss-Seidel //
  // ====================== //

  if( Comm().MyPID() == 0 ) cout << endl << "- Gauss-Seidel (sym)" << endl;

  for( double omega=0.25 ; omega<1.5 ; omega+=0.25)
  {

    Time.ResetStartTime();

    Teuchos::ParameterList NewList(List_);
    ML_Epetra::SetDefaults(Defaults, NewList);
    NewList.set("output", 0);

    NewList.set("smoother: type", "symmetric Gauss-Seidel");
    NewList.set("smoother: damping", omega);
    NewList.set("smoother: sweeps", sweeps);
    yo = new ML_Epetra::MultiLevelPreconditioner(*RowMatrix_,NewList, true);
    assert( yo != 0 );
    solver.SetPrecOperator(yo);

    LHS.PutScalar(0.0);
    RHS.Random();

    solver.Iterate(MaxIters,Tol);
    solver.GetAllAztecStatus(status);
    sprintf(smoother,"n=%d, omega=%5.2e", sweeps, omega);
    ReqTime = Time.ElapsedTime();
    if( ReqTime < BestTime ) BestTimeCount = count;
    if( (int) status[AZ_its] < BestIters ) BestItersCount = count;
    if( Comm().MyPID() == 0 ) MLP_print(count++,smoother,status,ReqTime);

    delete yo;
  }

  // ================== //
  // block Gauss-Seidel //
  // =================== //

  if( Comm().MyPID() == 0 ) cout << endl << "- Gauss-Seidel (block)" << endl;

  for( double omega=0.25 ; omega<1.5 ; omega+=0.25)
  {
    Time.ResetStartTime();

    Teuchos::ParameterList NewList(List_);
    ML_Epetra::SetDefaults(Defaults, NewList);
    NewList.set("output", 0);

    NewList.set("smoother: type", "block Gauss-Seidel");
    yo = new ML_Epetra::MultiLevelPreconditioner(*RowMatrix_,NewList, true);
    assert( yo != 0 );
    solver.SetPrecOperator(yo);

    LHS.PutScalar(0.0);
    RHS.Random();

    solver.Iterate(MaxIters,Tol);
    solver.GetAllAztecStatus(status);
    sprintf(smoother,"n=%d, omega=%5.2e", sweeps, omega);
    ReqTime = Time.ElapsedTime();
    if( ReqTime < BestTime ) BestTimeCount = count;
    if( (int) status[AZ_its] < BestIters ) BestItersCount = count;
    if( Comm().MyPID() == 0 ) MLP_print(count++,smoother,status,ReqTime);

    delete yo;
  } 

  // ==================== //
  // Aztec preconditioner //
  // ==================== //
  
  if( Comm().MyPID() == 0 ) cout << endl << "- Aztec preconditioner" << endl;
  
  for( int fillin=0 ; fillin<3 ; ++fillin ) {

    int options[AZ_OPTIONS_SIZE];
    double params[AZ_PARAMS_SIZE];
    AZ_defaults(options,params);
    options[AZ_graph_fill] = fillin;
    options[AZ_precond] = AZ_dom_decomp;
    options[AZ_subdomain_solve] = AZ_ilu;
    Time.ResetStartTime();

    Teuchos::ParameterList NewList(List_);
    ML_Epetra::SetDefaults(Defaults, NewList);
    NewList.set("output", 0);

    NewList.set("smoother: type", "Aztec");
    NewList.set("smoother: sweeps", sweeps);
    NewList.set("smoother: Aztec options", options);
    NewList.set("smoother: Aztec params", params);

    yo = new ML_Epetra::MultiLevelPreconditioner(*RowMatrix_,NewList, true);
    assert( yo != 0 );
    solver.SetPrecOperator(yo);

    LHS.PutScalar(0.0);
    RHS.Random();

    solver.Iterate(MaxIters,Tol);
    solver.GetAllAztecStatus(status);
    sprintf(smoother,"ILU(fill=%d)",fillin);
    ReqTime = Time.ElapsedTime();
    if( ReqTime < BestTime ) BestTimeCount = count;
    if( (int) status[AZ_its] < BestIters ) BestItersCount = count;
    if( Comm().MyPID() == 0 ) MLP_print(count++,smoother,status,ReqTime);

    delete yo;
  }

  // ================= //
  // Aztec as a solver //
  // ================= //
  
  if( Comm().MyPID() == 0 ) cout << endl << "- Aztec as solver" << endl;

  for( int iters=1 ; iters<6 ; iters+=2 ) {
    Time.ResetStartTime();

    Teuchos::ParameterList NewList(List_);
    ML_Epetra::SetDefaults(Defaults, NewList);
    NewList.set("output", 0);

    NewList.set("smoother: type", "Aztec");
    NewList.set("smoother: sweeps", iters);
    NewList.set("smoother: Aztec as solver", true);
    yo = new ML_Epetra::MultiLevelPreconditioner(*RowMatrix_,NewList, true);
    assert( yo != 0 );
     
    solver.SetAztecOption(AZ_solver, AZ_GMRESR); 
    solver.SetPrecOperator(yo);

    LHS.PutScalar(0.0);
    RHS.Random();

    solver.Iterate(MaxIters,Tol);
    solver.GetAllAztecStatus(status);
    sprintf(smoother,"iterations=%d", iters);
    ReqTime = Time.ElapsedTime();
    if( ReqTime < BestTime ) BestTimeCount = count;
    if( (int) status[AZ_its] < BestIters ) BestItersCount = count;
    if( Comm().MyPID() == 0 ) MLP_print(count++,smoother,status,ReqTime);
    solver.SetAztecOption(AZ_solver, AZ_gmres); 

    delete yo;
  }

  // ========= //
  // ParaSails //
  // ========= //
#ifdef HAVE_ML_PARASAILS
  if( Comm().MyPID() == 0 ) cout << endl << "- ParaSails" << endl;

  {
    
    Time.ResetStartTime();

    Teuchos::ParameterList NewList(List_);
    ML_Epetra::SetDefaults(Defaults, NewList);
    NewList.set("output", 0);

    NewList.set("smoother: type", "ParaSails");
    yo = new ML_Epetra::MultiLevelPreconditioner(*RowMatrix_,NewList, true);
    assert( yo != 0 );
    solver.SetPrecOperator(yo);

    LHS.PutScalar(0.0);
    RHS.Random();

    solver.Iterate(MaxIters,Tol);
    solver.GetAllAztecStatus(status);
    sprintf(smoother,"default");
    ReqTime = Time.ElapsedTime();
    if( ReqTime < BestTime ) BestTimeCount = count;
    if( (int) status[AZ_its] < BestIters ) BestItersCount = count;
    if( Comm().MyPID() == 0 ) MLP_print(count++,smoother,status,ReqTime);

    delete yo;
  }
#endif

  // ====== //
  // IFPACK //
  // ====== //

#ifdef HAVE_ML_IFPACKzzz
  if( Comm().MyPID() == 0 ) cout << endl << "- IFPACK" << endl;

  {
    Time.ResetStartTime();

    Teuchos::ParameterList NewList(List_);
    ML_Epetra::SetDefaults(Defaults, NewList);
    NewList.set("output", 0);

    NewList.set("smoother: type", "IFPACK");
    yo = new ML_Epetra::MultiLevelPreconditioner(*RowMatrix_,NewList, true);
    assert( yo != 0 );
     
    solver.SetPrecOperator(yo);

    LHS.PutScalar(0.0);
    RHS.Random();

    solver.Iterate(MaxIters,Tol);
    solver.GetAllAztecStatus(status);
    sprintf(smoother,"default");
    ReqTime = Time.ElapsedTime();
    if( ReqTime < BestTime ) BestTimeCount = count;
    if( (int) status[AZ_its] < BestIters ) BestItersCount = count;
    if( Comm().MyPID() == 0 ) MLP_print(count++,smoother,status,ReqTime);

    delete yo;
  }
#endif

  if( Comm().MyPID() == 0 ) {
    cout << "@@@ The best iteration count was obtain in test " << BestItersCount << endl;
    cout << "@@@ The best CPU-time was obtain in test " << BestTimeCount << endl;
  }

  // ================ //
  // that's all folks //
  // ================ //
  return;
}
// ============================================================================
void ML_Epetra::MultiLevelPreconditioner::AnalyzeMatrixCoarsening(char * Defaults,
								  bool IsSymmetric)
{
  return;
}

#endif /*ifdef ML_WITH_EPETRA && ML_HAVE_TEUCHOS*/
