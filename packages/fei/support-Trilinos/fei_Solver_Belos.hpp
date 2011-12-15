#ifndef _fei_Solver_Belos_hpp_
#define _fei_Solver_Belos_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_trilinos_macros.hpp>

#ifdef HAVE_FEI_BELOS

#include <fei_macros.hpp>
#include <fei_Solver.hpp>
#include <fei_Logger.hpp>

#ifdef HAVE_FEI_TEUCHOS
namespace Teuchos {
  class ParameterList;
}
#endif

class Epetra_CrsMatrix;
namespace Belos {
template<class Scalar,class MV,class OP> class SolverManager;
}

#ifdef HAVE_FEI_ML
#include <ml_include.h>
#include <ml_epetra_preconditioner.h>
#endif

/** fei::Solver implementation that wraps Trilinos/Belos.
 */
class Solver_Belos : public fei::Solver, private fei::Logger {
 public:
  /** Constructor
      No arguments required at construction time.
  */
  Solver_Belos();

  /** Destructor
   */
  virtual ~Solver_Belos();

  /** Method that accepts a fei::LinearSystem, extracts
      Epetra matrix/vector objects (A, x, b), and solves
      the linear-system using Belos. (After first attempting
      to parse control parameters that may specify things
      like tolerance, solver algorithm, etc.)

      @param numParams Number of solver-parameters (number
      of strings in the 'solverParams' list-of-strings)

      @param solverParams List of strings which are assumed
      to contain space-separated control parameters like
      "AZ_solver AZ_gmres", "AZ_tol 1.e-6", etc.
  */
  int solve(fei::LinearSystem* linearSystem,
	    fei::Matrix* preconditioningMatrix,
	    int numParams,
	    const char* const* solverParams,
	    int& iterationsTaken,
	    int& status);

  /** Method that accepts a fei::LinearSystem, extracts
      Epetra matrix/vector objects (A, x, b), and solves
      the linear-system using Belos. (After
      first attempting to parse control parameters that
      may specify things like tolerance, solver algorithm,
      etc.)

      @param parameterSet Object containing control
      parameters. It may have been populated using
      statements like:
      parameterSet.add(fei::Param("AZ_tol", 1.e-8));

      Note that if AZ_tol and/or AZ_max_iter are specified here,
      these values will override values that may have been
      set via the setMaxIters() and setTolerance() methods.
  */
  int solve(fei::LinearSystem* linearSystem,
	    fei::Matrix* preconditioningMatrix,
	    const fei::ParameterSet& parameterSet,
	    int& iterationsTaken,
	    int& status);

  Teuchos::ParameterList& get_ParameterList();

  void setMaxIters(int maxits) {maxIters_ = maxits;}
  void setTolerance(double tol) {tolerance_ = tol;}
  void setUseTranspose(bool useTrans) {useTranspose_ = useTrans;}

  void setUseML(bool useml);

 private:
  int setup_ml_operator(AztecOO& azoo, Epetra_CrsMatrix* A);

 private:
  Teuchos::RCP<Belos::SolverManager<double,Epetra_MultiVector,Epetra_Operator> > belos_solver_manager_;
  double tolerance_;
  int maxIters_;
  bool useTranspose_;
  Teuchos::RCP<Teuchos::ParameterList> paramlist_;

  bool useML_;
#ifdef HAVE_FEI_ML
  ML_Epetra::MultiLevelPreconditioner* ml_prec_;
  bool ml_defaults_set_;
  int *ml_aztec_options_;
  double *ml_aztec_params_;
#endif

  std::string name_;
  std::string dbgprefix_;
}; //class Solver_Belos

#endif // HAVE_FEI_BELOS

#endif
