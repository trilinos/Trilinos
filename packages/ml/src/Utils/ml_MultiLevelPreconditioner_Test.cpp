/*!
 *  \file ml_MultiLevelPreconditioner_Analyze.cpp
 *
 *  \brief Visualization utilities for MultiLevelPreconditioner class
 *
 *  \author Marzio Sala, SNL, 9214
 *
 *  \date Last update to Doxygen: 09-Aug-04
 *
 */
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

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
#include "ml_amesos_wrap.h"
#include "ml_MultiLevelPreconditioner.h"
#include "AztecOO.h"

// ============================================================================
static int SetLHSAndRHS(Epetra_Vector& LHS, Epetra_Vector& RHS,
			const Epetra_RowMatrix& Matrix)
{
  ML_CHK_ERR(LHS.Random());

  ML_CHK_ERR(Matrix.Multiply(false,LHS,RHS));

  ML_CHK_ERR(LHS.PutScalar(0.0));

  return(0);

}
// ============================================================================
static void MLP_print(int count, char * str, double status[AZ_STATUS_SIZE], double time)
{

  std::cout << "#" << count;
  if( count < 10  )      std::cout << ".....";
  else if( count < 100 ) std::cout << "....";
  else std::cout << "...";

  std::cout.width(30); std::cout.setf(std::ios::left);
  std::cout << str;
  std::cout.width(10); std::cout.setf(std::ios::left);
  std::cout << (int)status[AZ_its];
  std::cout.width(15); std::cout.setf(std::ios::left);
  std::cout << status[AZ_scaled_r];
  std::cout.width(15); std::cout.setf(std::ios::left);
  std::cout << time;
  
  if( status[AZ_why] == AZ_normal         ) std::cout << "N";
  else if( status[AZ_why] == AZ_maxits    ) std::cout << "M";
  else if( status[AZ_why] == AZ_loss      ) std::cout << "L";
  else if( status[AZ_why] == AZ_ill_cond  ) std::cout << "I";
  else if( status[AZ_why] == AZ_breakdown ) std::cout << "B";

  std::cout << std::endl;
}

// ============================================================================
int ML_Epetra::MultiLevelPreconditioner::
TestSmoothers(Teuchos::ParameterList& InputList,
	      const bool IsSymmetric)
{

  // sanity checks

  if (RowMatrix_ == 0) 
    ML_CHK_ERR(-1); // Matrix not yet set

  if (ML_isKLUAvailable() == 0)
    ML_CHK_ERR(-1); // need Amesos/KLU installed to do this

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
    std::cout << std::endl;
    ML_print_line("-",78);
    std::cout << "*** Analysis of ML parameters (smoothers)" << std::endl;
    std::cout << std::endl;;
    std::cout << "maximum iterations = " << MaxIters << std::endl;
    std::cout << "tolerance          = " << Tol << std::endl << std::endl;
    std::cout << "All options as in the input parameter list, except that" << std::endl;
    std::cout << "all levels have the same smoother" << std::endl << std::endl;
    std::cout << "M: maximum iterations exceeded without convergence" << std::endl;
    std::cout << "N: normal exit status (convergence achieved)" << std::endl;
    std::cout << "B: breakdown occurred" << std::endl;
    std::cout << "I: matrix is ill-conditioned" << std::endl;
    std::cout << "L: numerical loss of precision occurred" << std::endl;
    std::cout << std::endl;
    std::cout << "count  ";
    std::cout.width(30); std::cout.setf(std::ios::left); 
    std::cout << "smoother type";
    std::cout.width(10); std::cout.setf(std::ios::left);
    std::cout << "its";
    std::cout.width(15); std::cout.setf(std::ios::left);
    std::cout << "||r||/||r_0||";
    std::cout.width(15); std::cout.setf(std::ios::left);
    std::cout << "time (s)" << std::endl;
  }

  // FIXME: this use of LevelID_ is dangerous...
  // better to replace with something else
  // ====== //
  // Jacobi //
  // ====== //

  if (InputList.get("test: Jacobi",true) == true) {

    if( Comm().MyPID() == 0 ) std::cout << std::endl << "- Jacobi" << std::endl;

    for( double omega=0.25 ; omega<1.5 ; omega+=0.25) {

      Time.ResetStartTime();

      Teuchos::ParameterList NewList(InputList);
      NewList.set("ML output", 0);
      for (int ilevel = 0 ; ilevel < MaxLevels ; ++ilevel) {
	sprintf(parameter,"smoother: type (level %d)", LevelID_[ilevel]);
	NewList.set(parameter, "Jacobi");
	sprintf(parameter,"smoother: damping factor (level %d)", LevelID_[ilevel]);
	NewList.set(parameter, omega);
	sprintf(parameter,"smoother: sweeps (level %d)", LevelID_[ilevel]);
	NewList.set(parameter, sweeps);
      }
      yo = new ML_Epetra::MultiLevelPreconditioner(*RowMatrix_,NewList, true);
      assert( yo != 0 );
      solver.SetPrecOperator(yo);

      SetLHSAndRHS(LHS, RHS, *RowMatrix_);

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

    if( Comm().MyPID() == 0 ) std::cout << std::endl << "- Gauss-Seidel" << std::endl;

    for( double omega=0.25 ; omega<1.5 ; omega+=0.25) {

      Time.ResetStartTime();

      Teuchos::ParameterList NewList(InputList);
      NewList.set("ML output", 0);

      for (int ilevel = 0 ; ilevel < MaxLevels ; ++ilevel) {
	sprintf(parameter,"smoother: type (level %d)", LevelID_[ilevel]);
	NewList.set(parameter, "Gauss-Seidel");
	sprintf(parameter,"smoother: damping factor (level %d)", LevelID_[ilevel]);
	NewList.set(parameter, omega);
	sprintf(parameter,"smoother: sweeps (level %d)", LevelID_[ilevel]);
	NewList.set(parameter, sweeps);
      }
      yo = new ML_Epetra::MultiLevelPreconditioner(*RowMatrix_,NewList, true);
      assert( yo != 0 );
      solver.SetPrecOperator(yo);

      SetLHSAndRHS(LHS, RHS, *RowMatrix_);

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

    if( Comm().MyPID() == 0 ) std::cout << std::endl << "- Gauss-Seidel (sym)" << std::endl;

    for( double omega=0.25 ; omega<1.5 ; omega+=0.25) {

      Time.ResetStartTime();

      Teuchos::ParameterList NewList(InputList);
      NewList.set("ML output", 0);

      for (int ilevel = 0 ; ilevel < MaxLevels ; ++ilevel) {
	sprintf(parameter,"smoother: type (level %d)", LevelID_[ilevel]);
	NewList.set(parameter, "symmetric Gauss-Seidel");
	sprintf(parameter,"smoother: damping factor (level %d)", LevelID_[ilevel]);
	NewList.set(parameter, omega);
	sprintf(parameter,"smoother: sweeps (level %d)", LevelID_[ilevel]);
	NewList.set(parameter, sweeps);
      }
      yo = new ML_Epetra::MultiLevelPreconditioner(*RowMatrix_,NewList, true);
      assert( yo != 0 );
      solver.SetPrecOperator(yo);

      SetLHSAndRHS(LHS, RHS, *RowMatrix_);

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

    if( Comm().MyPID() == 0 ) std::cout << std::endl << "- Gauss-Seidel (block)" << std::endl;

    for( double omega=0.25 ; omega<1.5 ; omega+=0.25)
    {
      Time.ResetStartTime();

      Teuchos::ParameterList NewList(InputList);
      NewList.set("ML output", 0);

      for (int ilevel = 0 ; ilevel < MaxLevels ; ++ilevel) {
	sprintf(parameter,"smoother: type (level %d)", LevelID_[ilevel]);
	NewList.set(parameter, "block Gauss-Seidel");
	sprintf(parameter,"smoother: damping factor (level %d)", LevelID_[ilevel]);
	NewList.set(parameter, omega);
	sprintf(parameter,"smoother: sweeps (level %d)", LevelID_[ilevel]);
	NewList.set(parameter, sweeps);
      }
      yo = new ML_Epetra::MultiLevelPreconditioner(*RowMatrix_,NewList, true);
      assert( yo != 0 );
      solver.SetPrecOperator(yo);

      SetLHSAndRHS(LHS, RHS, *RowMatrix_);

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

    if( Comm().MyPID() == 0 ) std::cout << std::endl << "- Aztec preconditioner" << std::endl;

    for( int fillin=0 ; fillin<3 ; ++fillin ) {

      Teuchos::RCP<std::vector<int> > aztecOptions = Teuchos::rcp(new std::vector<int>(AZ_OPTIONS_SIZE));
      Teuchos::RCP<std::vector<double> > aztecParams  = Teuchos::rcp(new std::vector<double>(AZ_PARAMS_SIZE));

      int* options=  &(*aztecOptions)[0];
      double* params =  &(*aztecParams)[0];

      AZ_defaults(options,params);
      options[AZ_graph_fill] = fillin;
      options[AZ_precond] = AZ_dom_decomp;
      options[AZ_subdomain_solve] = AZ_ilu;
      Time.ResetStartTime();

      Teuchos::ParameterList NewList(InputList);
      NewList.set("ML output", 0);

      for (int ilevel = 0 ; ilevel < MaxLevels ; ++ilevel) {
	sprintf(parameter,"smoother: type (level %d)", LevelID_[ilevel]);
	NewList.set(parameter, "Aztec");
	sprintf(parameter,"smoother: sweeps (level %d)", LevelID_[ilevel]);
	NewList.set(parameter, sweeps);
	sprintf(parameter,"smoother: Aztec options (level %d)", LevelID_[ilevel]);
	NewList.set(parameter, aztecOptions);
	sprintf(parameter,"smoother: Aztec params (level %d)", LevelID_[ilevel]);
	NewList.set(parameter, aztecParams);
      }

      yo = new ML_Epetra::MultiLevelPreconditioner(*RowMatrix_,NewList, true);
      assert( yo != 0 );
      solver.SetPrecOperator(yo);

      SetLHSAndRHS(LHS, RHS, *RowMatrix_);

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

    if( Comm().MyPID() == 0 ) std::cout << std::endl << "- Aztec as solver" << std::endl;

    for( int iters=1 ; iters<6 ; iters+=2 ) {
      Time.ResetStartTime();

      Teuchos::ParameterList NewList(InputList);
      NewList.set("ML output", 0);

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

      SetLHSAndRHS(LHS, RHS, *RowMatrix_);

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

    if( Comm().MyPID() == 0 ) std::cout << std::endl << "- ParaSails" << std::endl;

    Time.ResetStartTime();

    Teuchos::ParameterList NewList(InputList);
    NewList.set("ML output", 0);

    for (int ilevel = 0 ; ilevel < MaxLevels ; ++ilevel) {
      sprintf(parameter,"smoother: type (level %d)", LevelID_[ilevel]);
      NewList.set(parameter, "ParaSails");
    }
    yo = new ML_Epetra::MultiLevelPreconditioner(*RowMatrix_,NewList, true);
    assert( yo != 0 );
    solver.SetPrecOperator(yo);

    SetLHSAndRHS(LHS, RHS, *RowMatrix_);

    solver.Iterate(MaxIters,Tol);
    solver.GetAllAztecStatus(status);
    sprintf(smoother,"%s","default");
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

#if defined(HAVE_ML_IFPACK) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)
  if (InputList.get("test: IFPACK",true) == true) {
    
    if( Comm().MyPID() == 0 ) std::cout << std::endl << "- IFPACK" << std::endl;

    // test IFPACK with block symmetric Gauss-Seidel, block
    // defined by the aggregates, and dense.

    for (double omega = 0.25 ; omega < 1.5 ; omega += 0.25) {

      Time.ResetStartTime();

      Teuchos::ParameterList NewList(InputList);
      NewList.set("ML output", 0);

      for (int ilevel = 0 ; ilevel < MaxLevels ; ++ilevel) {
        sprintf(parameter,"smoother: type (level %d)", LevelID_[ilevel]);
        NewList.set(parameter, "IFPACK");
      }

      NewList.set("smoother: ifpack type", "block relaxation stand-alone");
      NewList.set("smoother: ifpack overlap", 0);

      Teuchos::ParameterList& IFPACKList = NewList.sublist("smoother: ifpack list");
      IFPACKList.set("partitioner: type", "user");
      IFPACKList.set("relaxation: type", "symmetric Gauss-Seidel");
      IFPACKList.set("relaxation: damping factor", omega);
      IFPACKList.set("relaxation: sweeps", 1);
      IFPACKList.set("relaxation: zero starting solution", false);

      yo = new ML_Epetra::MultiLevelPreconditioner(*RowMatrix_,NewList, true);
      assert( yo != 0 );

      solver.SetPrecOperator(yo);

      SetLHSAndRHS(LHS, RHS, *RowMatrix_);

      solver.Iterate(MaxIters,Tol);
      solver.GetAllAztecStatus(status);
      sprintf(smoother,"BSGS-a, omega=%5.2e", omega);
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

#if defined(HAVE_ML_AMESOS)

    for (double omega = 0.25 ; omega < 1.5 ; omega += 0.25) {

      // now test block w/ equation partitioner

      Time.ResetStartTime();

      Teuchos::ParameterList NewList(InputList);
      NewList.set("ML output", 0);

      for (int ilevel = 0 ; ilevel < MaxLevels ; ++ilevel) {
        sprintf(parameter,"smoother: type (level %d)", LevelID_[ilevel]);
        NewList.set(parameter, "IFPACK");
      }

      NewList.set("smoother: ifpack type", "block relaxation stand-alone (Amesos)");
      NewList.set("smoother: ifpack overlap", 0);

      Teuchos::ParameterList& IFPACKList = NewList.sublist("smoother: ifpack list");
      IFPACKList.set("partitioner: type", "equation");
      IFPACKList.set("partitioner: local parts", NumPDEEqns_);
      IFPACKList.set("relaxation: type", "symmetric Gauss-Seidel");
      IFPACKList.set("relaxation: damping factor", omega);
      IFPACKList.set("relaxation: sweeps", 1);
      IFPACKList.set("relaxation: zero starting solution", false);

      yo = new ML_Epetra::MultiLevelPreconditioner(*RowMatrix_,NewList, true);
      assert( yo != 0 );

      solver.SetPrecOperator(yo);

      SetLHSAndRHS(LHS, RHS, *RowMatrix_);

      solver.Iterate(MaxIters,Tol);
      solver.GetAllAztecStatus(status);
      sprintf(smoother,"BSGS-e, omega=%5.2e", omega);
      ReqTime = Time.ElapsedTime();
      if (ReqTime < BestTime) {
        BestTime = ReqTime;
        BestTimeCount = count;
      }
      if ((int) status[AZ_its] < BestIters) {
        BestIters = (int)status[AZ_its];
        BestItersCount = count;
      }
      if (Comm().MyPID() == 0) MLP_print(count++,smoother,status,ReqTime);

      delete yo;
    }
#endif
 }
#endif

  // ================ //
  // ML self smoother //
  // ================ //

  if (InputList.get("test: ML self smoother",true) == true) {

    if( Comm().MyPID() == 0 ) std::cout << std::endl << "- ML as local smoother" << std::endl;

    Time.ResetStartTime();

    Teuchos::ParameterList NewList(InputList);
    NewList.set("ML output", 0);

    for (int ilevel = 0 ; ilevel < MaxLevels ; ++ilevel) {
      sprintf(parameter,"smoother: type (level %d)", LevelID_[ilevel]);
      NewList.set(parameter, "self");
    }

    Teuchos::ParameterList& SelfList = NewList.sublist("smoother: self list");
    SetDefaults("DD-ML", SelfList);
    SelfList.set("ML output", 0);

    yo = new ML_Epetra::MultiLevelPreconditioner(*RowMatrix_,NewList, true);
    assert (yo != 0);

    solver.SetPrecOperator(yo);

    SetLHSAndRHS(LHS, RHS, *RowMatrix_);

    solver.Iterate(MaxIters,Tol);
    solver.GetAllAztecStatus(status);
    sprintf(smoother,"%s","ML self, DD-ML");
    ReqTime = Time.ElapsedTime();
    if (ReqTime < BestTime) {
      BestTime = ReqTime;
      BestTimeCount = count;
    }
    if ((int) status[AZ_its] < BestIters) {
      BestIters = (int)status[AZ_its];
      BestItersCount = count;
    }
    if (Comm().MyPID() == 0) MLP_print(count++,smoother,status,ReqTime);

    delete yo;

  }

  if (Comm().MyPID() == 0) {
    std::cout << std::endl;
    std::cout << "*** The best iteration count was obtain in test " << BestItersCount << std::endl;
    std::cout << "*** The best CPU-time was obtain in test " << BestTimeCount << std::endl;
    std::cout << std::endl;
  }

  // ================ //
  // that's all folks //
  // ================ //

  if (Comm().MyPID() == 0) {
    std::cout << std::endl << "*** Total time = " << GetClock() - time << "(s)" << std::endl;
    ML_print_line("-",78);
    std::cout << std::endl;
  }

  return(0);
}

#endif /*ifdef ML_WITH_EPETRA && ML_HAVE_TEUCHOS*/
