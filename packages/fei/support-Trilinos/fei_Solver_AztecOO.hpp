#ifndef _fei_Solver_AztecOO_h_
#define _fei_Solver_AztecOO_h_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifdef HAVE_FEI_AZTECOO

#include <fei_macros.hpp>
#include <fei_Solver.hpp>
#include <fei_Logger.hpp>

#ifdef HAVE_FEI_TEUCHOS
namespace Teuchos {
  class ParameterList;
}
#endif

class AztecOO;

#ifdef HAVE_FEI_ML
#include <ml_include.h>
#include <ml_epetra_preconditioner.h>
#endif

/** fei::Solver implementation that wraps Trilinos/AztecOO.
 */
class Solver_AztecOO : public fei::Solver, private fei::Logger {
 public:
  /** Constructor
      No arguments required at construction time.
  */
  Solver_AztecOO();

  /** Destructor
   */
  virtual ~Solver_AztecOO();

  /** Method that accepts a fei::LinearSystem, extracts
      Epetra matrix/vector objects (A, x, b), and solves
      the linear-system using AztecOO::Iterate. (After
      first attempting to parse control parameters that
      may specify things like tolerance, solver algorithm,
      etc.)

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
      the linear-system using AztecOO::Iterate. (After
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

  AztecOO& getAztecOO();

  Teuchos::ParameterList& get_ParameterList();

  void setMaxIters(int maxits) {maxIters_ = maxits;}
  void setTolerance(double tol) {tolerance_ = tol;}
  void setUseTranspose(bool useTrans) {useTranspose_ = useTrans;}

  void setUseML(bool useml);

 private:
  int setup_ml_operator(AztecOO& azoo, Epetra_CrsMatrix* A);

 private:
  AztecOO* azoo_;
  double tolerance_;
  int maxIters_;
  bool useTranspose_;
  Teuchos::ParameterList* paramlist_;

  bool useML_;
#ifdef HAVE_FEI_ML
  ML_Epetra::MultiLevelPreconditioner* ml_prec_;
  bool ml_defaults_set_;
  int *ml_aztec_options_;
  double *ml_aztec_params_;
#endif

  std::string name_;
  std::string dbgprefix_;
}; //class Solver_AztecOO

#endif // HAVE_FEI_AZTECOO

#endif
