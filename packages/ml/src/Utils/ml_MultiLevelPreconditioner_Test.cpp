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
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_AZTECOO)

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Time.h"
#include "Epetra_RowMatrix.h"
#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "ml_epetra.h"
#include "ml_MultiLevelPreconditioner.h"
#include "AztecOO.h"

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
int ML_Epetra::MultiLevelPreconditioner::
TestSmoothers(Teuchos::ParameterList& InputList,
	      const bool IsSymmetric)
{
 
  // sanity checks

  if (RowMatrix_ == 0) 
    ML_CHK_ERR(-1); // Matrix not yet set

  // execution begins

  double time;

  time = GetClock();

  int MaxIters = InputList.get("test: max iters",500);
  double Tol   = InputList.get("test: tolerance", 1e-5);
  double status[AZ_STATUS_SIZE];
  char smoother[80];
  int sweeps = InputList.get("test: sweeps",1);
  Epetra_Time Time(Comm());
  int count = 0;
  double ReqTime;
  int BestIters = 1000000;
  double BestTime = 1000000.0;
  int BestItersCount = -1, BestTimeCount = -1;
  char parameter[80];
  int MaxLevels = InputList.get("max levels",2);

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
    cout << endl;
    ML_print_line("-",78);
    cout << "*** ************************************* ***" << endl;
    cout << "*** Analysis of ML parameters (smoothers) ***" << endl;
    cout << "*** ************************************* ***" << endl;
    cout << endl;;
    cout << "*** maximum iterations = " << MaxIters << endl;
    cout << "*** tolerance          = " << Tol << endl << endl;
    cout << "*** All options as in the input parameter list, except that" << endl;
    cout << "*** all levels have the same smoother" << endl << endl;
    cout << "*** M: maximum iterations exceeded without convergence" << endl;
    cout << "*** N: normal exit status (convergence achieved)" << endl;
    cout << "*** B: breakdown occurred" << endl;
    cout << "*** I: matrix is ill-conditioned" << endl;
    cout << "*** L: numerical loss of precision occurred" << endl;
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

  // FIXME: this use of LevelID_ is dangerous...
  // better to replace with something else
  // ====== //
  // Jacobi //
  // ====== //

  if (InputList.get("test: Jacobi",true) == true) {

    if( Comm().MyPID() == 0 ) cout << endl << "- Jacobi" << endl;

    for( double omega=0.25 ; omega<1.5 ; omega+=0.25) {

      Time.ResetStartTime();

      Teuchos::ParameterList NewList(InputList);
      NewList.set("output", 0);
      for (int ilevel = 0 ; ilevel < MaxLevels ; ++ilevel) {
	sprintf(parameter,"smoother: type (level %d)", LevelID_[ilevel]);
	NewList.set(parameter, "Jacobi");
	sprintf(parameter,"smoother: damping (level %d)", LevelID_[ilevel]);
	NewList.set(parameter, omega);
	sprintf(parameter,"smoother: sweeps (level %d)", LevelID_[ilevel]);
	NewList.set(parameter, sweeps);
      }
      yo = new ML_Epetra::MultiLevelPreconditioner(*RowMatrix_,NewList, true);
      assert( yo != 0 );
      solver.SetPrecOperator(yo);

      LHS.PutScalar(0.0);
      RHS.Random();

      solver.Iterate(MaxIters,Tol);
      solver.GetAllAztecStatus(status);
      sprintf(smoother,"n=%d, omega=%5.2e", sweeps, omega);
      ReqTime = Time.ElapsedTime();
      if (ReqTime < BestTime) {
	BestTime = ReqTime;
	BestTimeCount = count;
      }
      if ((int) status[AZ_its] < BestIters) {
	BestIters = (int)status[AZ_its];
	BestItersCount = count;
      }
      if( Comm().MyPID() == 0 ) MLP_print(count++,smoother,status,ReqTime);

      delete yo;
    }

  }

  // ============ //
  // Gauss-Seidel //
  // ============ //

  if (InputList.get("test: Gauss-Seidel",true) == true) {

    if( Comm().MyPID() == 0 ) cout << endl << "- Gauss-Seidel" << endl;

    for( double omega=0.25 ; omega<1.5 ; omega+=0.25) {

      Time.ResetStartTime();

      Teuchos::ParameterList NewList(InputList);
      NewList.set("output", 0);

      for (int ilevel = 0 ; ilevel < MaxLevels ; ++ilevel) {
	sprintf(parameter,"smoother: type (level %d)", LevelID_[ilevel]);
	NewList.set(parameter, "Gauss-Seidel");
	sprintf(parameter,"smoother: damping (level %d)", LevelID_[ilevel]);
	NewList.set(parameter, omega);
	sprintf(parameter,"smoother: sweeps (level %d)", LevelID_[ilevel]);
	NewList.set(parameter, sweeps);
      }
      yo = new ML_Epetra::MultiLevelPreconditioner(*RowMatrix_,NewList, true);
      assert( yo != 0 );
      solver.SetPrecOperator(yo);

      LHS.PutScalar(0.0);
      RHS.Random();

      solver.Iterate(MaxIters,Tol);
      solver.GetAllAztecStatus(status);
      sprintf(smoother,"n=%d, omega=%5.2e", sweeps, omega);
      ReqTime = Time.ElapsedTime();
      if (ReqTime < BestTime) {
	BestTime = ReqTime;
	BestTimeCount = count;
      }
      if ((int) status[AZ_its] < BestIters) {
	BestIters = (int)status[AZ_its];
	BestItersCount = count;
      }
      if( Comm().MyPID() == 0 ) MLP_print(count++,smoother,status,ReqTime);

      delete yo;
    }
  }

  // ====================== //
  // symmetric Gauss-Seidel //
  // ====================== //

  if (InputList.get("test: symmetric Gauss-Seidel",true) == true) {

    if( Comm().MyPID() == 0 ) cout << endl << "- Gauss-Seidel (sym)" << endl;

    for( double omega=0.25 ; omega<1.5 ; omega+=0.25)
    {

      Time.ResetStartTime();

      Teuchos::ParameterList NewList(InputList);
      NewList.set("output", 0);

      for (int ilevel = 0 ; ilevel < MaxLevels ; ++ilevel) {
	sprintf(parameter,"smoother: type (level %d)", LevelID_[ilevel]);
	NewList.set(parameter, "symmetric Gauss-Seidel");
	sprintf(parameter,"smoother: damping (level %d)", LevelID_[ilevel]);
	NewList.set(parameter, omega);
	sprintf(parameter,"smoother: sweeps (level %d)", LevelID_[ilevel]);
	NewList.set(parameter, sweeps);
      }
      yo = new ML_Epetra::MultiLevelPreconditioner(*RowMatrix_,NewList, true);
      assert( yo != 0 );
      solver.SetPrecOperator(yo);

      LHS.PutScalar(0.0);
      RHS.Random();

      solver.Iterate(MaxIters,Tol);
      solver.GetAllAztecStatus(status);
      sprintf(smoother,"n=%d, omega=%5.2e", sweeps, omega);
      ReqTime = Time.ElapsedTime();
      if (ReqTime < BestTime) {
	BestTime = ReqTime;
	BestTimeCount = count;
      }
      if ((int) status[AZ_its] < BestIters) {
	BestIters = (int)status[AZ_its];
	BestItersCount = count;
      }
      if( Comm().MyPID() == 0 ) MLP_print(count++,smoother,status,ReqTime);

      delete yo;
    }
  }

  // ================== //
  // block Gauss-Seidel //
  // =================== //

  if (InputList.get("test: block Gauss-Seidel",true) == true) {

    if( Comm().MyPID() == 0 ) cout << endl << "- Gauss-Seidel (block)" << endl;

    for( double omega=0.25 ; omega<1.5 ; omega+=0.25)
    {
      Time.ResetStartTime();

      Teuchos::ParameterList NewList(InputList);
      NewList.set("output", 0);

      for (int ilevel = 0 ; ilevel < MaxLevels ; ++ilevel) {
	sprintf(parameter,"smoother: type (level %d)", LevelID_[ilevel]);
	NewList.set(parameter, "block Gauss-Seidel");
	sprintf(parameter,"smoother: damping (level %d)", LevelID_[ilevel]);
	NewList.set(parameter, omega);
	sprintf(parameter,"smoother: sweeps (level %d)", LevelID_[ilevel]);
	NewList.set(parameter, sweeps);
      }
      yo = new ML_Epetra::MultiLevelPreconditioner(*RowMatrix_,NewList, true);
      assert( yo != 0 );
      solver.SetPrecOperator(yo);

      LHS.PutScalar(0.0);
      RHS.Random();

      solver.Iterate(MaxIters,Tol);
      solver.GetAllAztecStatus(status);
      sprintf(smoother,"n=%d, omega=%5.2e", sweeps, omega);
      ReqTime = Time.ElapsedTime();
      if (ReqTime < BestTime) {
	BestTime = ReqTime;
	BestTimeCount = count;
      }
      if ((int) status[AZ_its] < BestIters) {
	BestIters = (int)status[AZ_its];
	BestItersCount = count;
      }
      if( Comm().MyPID() == 0 ) MLP_print(count++,smoother,status,ReqTime);

      delete yo;
    } 
  }

  // ==================== //
  // Aztec preconditioner //
  // ==================== //
  
  if (InputList.get("test: Aztec",true) == true) {

    if( Comm().MyPID() == 0 ) cout << endl << "- Aztec preconditioner" << endl;

    for( int fillin=0 ; fillin<3 ; ++fillin ) {

      int options[AZ_OPTIONS_SIZE];
      double params[AZ_PARAMS_SIZE];
      AZ_defaults(options,params);
      options[AZ_graph_fill] = fillin;
      options[AZ_precond] = AZ_dom_decomp;
      options[AZ_subdomain_solve] = AZ_ilu;
      Time.ResetStartTime();

      Teuchos::ParameterList NewList(InputList);
      NewList.set("output", 0);

      for (int ilevel = 0 ; ilevel < MaxLevels ; ++ilevel) {
	sprintf(parameter,"smoother: type (level %d)", LevelID_[ilevel]);
	NewList.set(parameter, "Aztec");
	sprintf(parameter,"smoother: sweeps (level %d)", LevelID_[ilevel]);
	NewList.set(parameter, sweeps);
	sprintf(parameter,"smoother: Aztec options (level %d)", LevelID_[ilevel]);
	NewList.set(parameter, options);
	sprintf(parameter,"smoother: Aztec params (level %d)", LevelID_[ilevel]);
	NewList.set(parameter, params);
      }

      yo = new ML_Epetra::MultiLevelPreconditioner(*RowMatrix_,NewList, true);
      assert( yo != 0 );
      solver.SetPrecOperator(yo);

      LHS.PutScalar(0.0);
      RHS.Random();

      solver.Iterate(MaxIters,Tol);
      solver.GetAllAztecStatus(status);
      sprintf(smoother,"ILU(fill=%d)",fillin);
      ReqTime = Time.ElapsedTime();
      if (ReqTime < BestTime) {
	BestTime = ReqTime;
	BestTimeCount = count;
      }
      if ((int) status[AZ_its] < BestIters) {
	BestIters = (int)status[AZ_its];
	BestItersCount = count;
      }
      if( Comm().MyPID() == 0 ) MLP_print(count++,smoother,status,ReqTime);

      delete yo;
    }
  }

  // ================= //
  // Aztec as a solver //
  // ================= //
  
  if (InputList.get("test: Aztec as solver",true) == true) {

    if( Comm().MyPID() == 0 ) cout << endl << "- Aztec as solver" << endl;

    for( int iters=1 ; iters<6 ; iters+=2 ) {
      Time.ResetStartTime();

      Teuchos::ParameterList NewList(InputList);
      NewList.set("output", 0);

      for (int ilevel = 0 ; ilevel < MaxLevels ; ++ilevel) {
	sprintf(parameter,"smoother: type (level %d)", LevelID_[ilevel]);
	NewList.set(parameter, "Aztec");
	sprintf(parameter,"smoother: Aztec as solver (level %d)", LevelID_[ilevel]);
	NewList.set(parameter, true);
	sprintf(parameter,"smoother: sweeps (level %d)", LevelID_[ilevel]);
	NewList.set(parameter, iters);
      }
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
      if (ReqTime < BestTime) {
	BestTime = ReqTime;
	BestTimeCount = count;
      }
      if ((int) status[AZ_its] < BestIters) {
	BestIters = (int)status[AZ_its];
	BestItersCount = count;
      }
      if( Comm().MyPID() == 0 ) MLP_print(count++,smoother,status,ReqTime);
      solver.SetAztecOption(AZ_solver, AZ_gmres); 

      delete yo;
    }
  }

  // ========= //
  // ParaSails //
  // ========= //
#ifdef HAVE_ML_PARASAILS
  if (InputList.get("test: ParaSails",true) == true) {

    if( Comm().MyPID() == 0 ) cout << endl << "- ParaSails" << endl;

    Time.ResetStartTime();

    Teuchos::ParameterList NewList(InputList);
    NewList.set("output", 0);

    for (int ilevel = 0 ; ilevel < MaxLevels ; ++ilevel) {
      sprintf(parameter,"smoother: type (level %d)", LevelID_[ilevel]);
      NewList.set(parameter, "ParaSails");
    }
    yo = new ML_Epetra::MultiLevelPreconditioner(*RowMatrix_,NewList, true);
    assert( yo != 0 );
    solver.SetPrecOperator(yo);

    LHS.PutScalar(0.0);
    RHS.Random();

    solver.Iterate(MaxIters,Tol);
    solver.GetAllAztecStatus(status);
    sprintf(smoother,"default");
    ReqTime = Time.ElapsedTime();
    if (ReqTime < BestTime) {
      BestTime = ReqTime;
      BestTimeCount = count;
    }
    if ((int) status[AZ_its] < BestIters) {
      BestIters = (int)status[AZ_its];
      BestItersCount = count;
    }
    if( Comm().MyPID() == 0 ) MLP_print(count++,smoother,status,ReqTime);

    delete yo;
  }
#endif

  // ====== //
  // IFPACK //
  // ====== //

  // FIXME: IFPACK is broken...
#ifdef HAVE_ML_IFPACKzzz
  if (InputList.get("test: IFPACK",true) == true) {
    
    if( Comm().MyPID() == 0 ) cout << endl << "- IFPACK" << endl;

    Time.ResetStartTime();

    Teuchos::ParameterList NewList(InputList);
    NewList.set("output", 0);

    for (int ilevel = 0 ; ilevel < MaxLevels ; ++ilevel) {
      sprintf(parameter,"smoother: type (level %d)", LevelID_[ilevel]);
      NewList.set(parameter, "IFPACK");
    }
    yo = new ML_Epetra::MultiLevelPreconditioner(*RowMatrix_,NewList, true);
    assert( yo != 0 );
     
    solver.SetPrecOperator(yo);

    LHS.PutScalar(0.0);
    RHS.Random();

    solver.Iterate(MaxIters,Tol);
    solver.GetAllAztecStatus(status);
    sprintf(smoother,"default");
    ReqTime = Time.ElapsedTime();
    if (ReqTime < BestTime) {
      BestTime = ReqTime;
      BestTimeCount = count;
    }
    if ((int) status[AZ_its] < BestIters) {
      BestIters = (int)status[AZ_its];
      BestItersCount = count;
    }
    if( Comm().MyPID() == 0 ) MLP_print(count++,smoother,status,ReqTime);

    delete yo;
  }
#endif

  if( Comm().MyPID() == 0 ) {
    cout << endl;
    cout << "*** The best iteration count was obtain in test " << BestItersCount << endl;
    cout << "*** The best CPU-time was obtain in test " << BestTimeCount << endl;
    cout << endl;
  }

  // ================ //
  // that's all folks //
  // ================ //

  if (Comm().MyPID() == 0) {
    cout << endl << "*** Total time = " << GetClock() - time << "(s)" << endl;
    ML_print_line("-",78);
    cout << endl;
  }

  return(0);
}

#endif /*ifdef ML_WITH_EPETRA && ML_HAVE_TEUCHOS*/
