#ifndef _fei_Solver_Amesos_h_
#define _fei_Solver_Amesos_h_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_trilinos_macros.hpp>

#ifdef HAVE_FEI_AMESOS

#include <fei_Solver.hpp>

namespace Teuchos {
  class ParameterList;
}
class Amesos;
class Amesos_BaseSolver;
class Epetra_LinearProblem;
class Epetra_CrsMatrix;
class Epetra_MultiVector;
class Epetra_VbrMatrix;

class Solver_Amesos : public fei::Solver {
 public:
  Solver_Amesos();
  virtual ~Solver_Amesos();

  int solve(fei::LinearSystem* linearSystem,
	    fei::Matrix* preconditioningMatrix,
	    const fei::ParameterSet& parameterSet,
	    int& iterationsTaken,
	    int& status);

  Teuchos::ParameterList& get_ParameterList();

 private:
  int solve(fei::LinearSystem* linearSystem,
	    fei::Matrix* preconditioningMatrix,
	    int numParams,
	    const char* const* solverParams,
	    int& iterationsTaken,
	    int& status);

  int parseParameters(int numParams,
		      const char*const* params);

  int solve_private(Epetra_CrsMatrix* A,
		    Epetra_MultiVector* x,
		    Epetra_MultiVector* b,
		    fei::Matrix* preconditioningMatrix,
		    int numParams,
		    const char* const* solverParams,
		    int& iterationsTaken,
		    int& status);

  int solve_private(Epetra_VbrMatrix* A,
		    Epetra_MultiVector* x,
		    Epetra_MultiVector* b,
		    fei::Matrix* preconditioningMatrix,
		    int numParams,
		    const char* const* solverParams,
		    int& iterationsTaken,
		    int& status);

 private:
  double tolerance_;
  int maxIters_;
  Amesos* amesos_factory_;
  Amesos_BaseSolver* amesos_solver_;
  Epetra_LinearProblem* epetra_linearproblem_;
  Teuchos::ParameterList* paramlist_;
}; //class Solver_Amesos

#endif

#endif
