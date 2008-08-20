/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <fei_Solver.hpp>

#include <fei_Matrix_Impl.hpp>
#include <fei_MatrixReducer.hpp>
#include <snl_fei_LinearSystem_FEData.hpp>
#include <fei_ParameterSet.hpp>
#include <fei_utils.hpp>

#undef fei_file
#define fei_file "fei_Solver.cpp"
#include <fei_ErrMacros.hpp>

//----------------------------------------------------------------------------
int fei_Solver_solve(fei::LinearSystem* linearSystem,
			   fei::Matrix* preconditioningMatrix,
			   int numParams,
			   const char* const* solverParams,
			   int& iterationsTaken,
			   int& status)
{
  fei::SharedPtr<fei::Matrix> matrix = linearSystem->getMatrix();
  fei::Matrix_Impl<LinearSystemCore>* lscmatrix =
    dynamic_cast<fei::Matrix_Impl<LinearSystemCore>*>(matrix.get());

  fei::MatrixReducer* matred = dynamic_cast<fei::MatrixReducer*>(matrix.get());
  if (matred != NULL) {
    lscmatrix = dynamic_cast<fei::Matrix_Impl<LinearSystemCore>*>(matred->getTargetMatrix().get());
  }

  if (lscmatrix != NULL) {
    fei::SharedPtr<LinearSystemCore> linSysCore = lscmatrix->getMatrix();

    char** params = const_cast<char**>(solverParams);
    CHK_ERR( linSysCore->parameters(numParams, params) );

    CHK_ERR( linSysCore->launchSolver(status, iterationsTaken) );

    return(0);
  }

  snl_fei::LinearSystem_FEData* fedlinsys =
    dynamic_cast<snl_fei::LinearSystem_FEData*>(linearSystem);
  if (fedlinsys != NULL) {
    fei::SharedPtr<FiniteElementData> fedata = fedlinsys->getFiniteElementData();

    CHK_ERR( fedata->launchSolver(status, iterationsTaken) );

    return(0);
  }

  ERReturn(-1);
}

//----------------------------------------------------------------------------
int fei::Solver::solve(fei::LinearSystem* linearSystem,
			   fei::Matrix* preconditioningMatrix,
			   const fei::ParameterSet& parameterSet,
			   int& iterationsTaken,
			   int& status)
{
  int numParams = 0;
  const char** paramStrings = NULL;
  std::vector<std::string> stdstrings;
  fei::utils::convert_ParameterSet_to_strings(&parameterSet, stdstrings);
  fei::utils::strings_to_char_ptrs(stdstrings, numParams, paramStrings);

  int err = fei_Solver_solve(linearSystem, preconditioningMatrix,
		  numParams, paramStrings,
		  iterationsTaken, status);

  delete [] paramStrings;

  return(err);
}

