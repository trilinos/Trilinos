/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_sstream.hpp>

#include "fei_trilinos_macros.hpp"

#ifdef HAVE_FEI_AMESOS
#ifdef HAVE_FEI_EPETRA

#include <fei_Solver_Amesos.hpp>
#include <fei_ParameterSet.hpp>
#include <snl_fei_Utils.hpp>
#include <fei_utils.hpp>

//fei_Include_Trilinos.hpp includes the actual Trilinos headers (epetra, aztecoo, ml ...)
#include <fei_Include_Trilinos.hpp>
#include <fei_Trilinos_Helpers.hpp>
#include <fei_VectorTraits_Epetra.hpp>
#include <fei_MatrixTraits_Epetra.hpp>

#include <fei_Vector.hpp>
#include <fei_Matrix.hpp>
#include <fei_LinearSystem.hpp>

//---------------------------------------------------------------------------
Solver_Amesos::Solver_Amesos()
  : tolerance_(1.e-6),
    maxIters_(500),
    amesos_factory_(new Amesos()),
    amesos_solver_(NULL),
    epetra_linearproblem_(NULL)
{
  paramlist_ = new Teuchos::ParameterList;
}

//---------------------------------------------------------------------------
Solver_Amesos::~Solver_Amesos()
{
  delete amesos_solver_;
  delete amesos_factory_;
  delete epetra_linearproblem_;
  delete paramlist_;
}

//---------------------------------------------------------------------------
Teuchos::ParameterList& Solver_Amesos::get_ParameterList()
{
  if (paramlist_ == NULL) {
    paramlist_ = new Teuchos::ParameterList;
  }

  return( *paramlist_ );
}

//---------------------------------------------------------------------------
int Solver_Amesos::solve(fei::LinearSystem* linearSystem,
			  fei::Matrix* preconditioningMatrix,
			  const fei::ParameterSet& parameterSet,
			  int& iterationsTaken,
			  int& status)
{
  Trilinos_Helpers::copy_parameterset(parameterSet, *paramlist_);

  int numParams = 0;
  const char** paramStrings = NULL;
  std::vector<std::string> stdstrings;
  fei::utils::convert_ParameterSet_to_strings(&parameterSet, stdstrings);
  fei::utils::strings_to_char_ptrs(stdstrings, numParams, paramStrings);

  int err = solve(linearSystem, preconditioningMatrix, numParams, paramStrings,
		  iterationsTaken, status);

  int olevel = 0;
  parameterSet.getIntParamValue("outputLevel", olevel);

  std::string param2;
  parameterSet.getStringParamValue("FEI_OUTPUT_LEVEL", param2);

  if (olevel >= 3 || param2 == "MATRIX_FILES") {
    std::string param1;
    parameterSet.getStringParamValue("debugOutput", param1);

    FEI_OSTRINGSTREAM osstr;
    if (!param1.empty()) {
      osstr << param1 << "/";
    }
    else osstr << "./";

    static int counter = 1;
    osstr << "x_Amesos.vec.slv"<<counter++;
    fei::SharedPtr<fei::Vector> feix = linearSystem->getSolutionVector();
    feix->writeToFile(osstr.str().c_str());
  }

  delete [] paramStrings;
  
  return(err);
}

//---------------------------------------------------------------------------
int Solver_Amesos::solve(fei::LinearSystem* linearSystem,
			  fei::Matrix* preconditioningMatrix,
			  int numParams,
			  const char* const* solverParams,
			  int& iterationsTaken,
			  int& status)
{
  Epetra_MultiVector*    x = NULL;
  Epetra_MultiVector*    b = NULL;
  Epetra_CrsMatrix* A = NULL;
  Epetra_Operator* opA = NULL;

  fei::SharedPtr<fei::Matrix> feiA = linearSystem->getMatrix();
  fei::SharedPtr<fei::Vector> feix = linearSystem->getSolutionVector();
  fei::SharedPtr<fei::Vector> feib = linearSystem->getRHS();

  Trilinos_Helpers::get_Epetra_pointers(feiA, feix, feib,
                                        A, opA, x, b);

  if (opA == 0 || x == 0 || b == 0) {
    fei::console_out() << "Error, couldn't obtain Epetra objects from "
      << "fei container-objects."<<FEI_ENDL;
    return(-1);
  }

  if (epetra_linearproblem_ == NULL) {
    epetra_linearproblem_ = new Epetra_LinearProblem;
  }

  epetra_linearproblem_->SetOperator(A);
  epetra_linearproblem_->SetLHS(x);
  epetra_linearproblem_->SetRHS(b);

  const char* param = snl_fei::getParamValue("Trilinos_Solver",
					     numParams, solverParams);
  if (param != NULL) {
    if (amesos_solver_ == 0) {
      amesos_solver_ = amesos_factory_->Create(param, *epetra_linearproblem_);
    }
    if (amesos_solver_ == 0) {
      cerr << "Solver_Amesos::solve ERROR, couldn't create Amesos solver named "
	   << param << ", amesos_factory::Create returned NULL." << endl;
      status = -999;
      return(-1);
    }
  }
  else {
    static char amesosklu[] = "Amesos_Klu";
    if (amesos_solver_ == 0) {
      amesos_solver_ = amesos_factory_->Create( amesosklu,
						*epetra_linearproblem_);
    }
    if (amesos_solver_ == 0) {
      cerr << "Solver_Amesos::solve ERROR, couldn't create Amesos solver named "
	   << amesosklu << ", it's apparently not supported." << endl;
      status = -999;
      return(-1);
    }
  }

  amesos_solver_->SetParameters(*paramlist_);

  if (feiA->changedSinceMark()) {
    amesos_solver_->SymbolicFactorization();
    amesos_solver_->NumericFactorization();
    feiA->markState();
  }

  amesos_solver_->Solve();
  status = 0;

  return(0);
}

//---------------------------------------------------------------------------
int Solver_Amesos::parseParameters(int numParams,
				   const char*const* params)
{
 
  return(0);
}

#endif //HAVE_FEI_EPETRA
#endif //HAVE_FEI_AMESOS
