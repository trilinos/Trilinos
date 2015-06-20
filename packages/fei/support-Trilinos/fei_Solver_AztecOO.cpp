/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


#include "fei_trilinos_macros.hpp"

#ifdef HAVE_FEI_AZTECOO

//fei_Include_Trilinos.h includes the actual Trilinos headers (epetra, aztecoo, ...)
#include <fei_Include_Trilinos.hpp>
#include <fei_Trilinos_Helpers.hpp>

#include <fei_Solver_AztecOO.hpp>
#include <fei_LinProbMgr_EpetraBasic.hpp>
#include <fei_ParameterSet.hpp>
#include <fei_utils.hpp>

#include <fei_VectorTraits_Epetra.hpp>
#include <fei_Vector_Impl.hpp>
#include <fei_MatrixTraits_Epetra.hpp>
#include <fei_Matrix_Impl.hpp>
#include <fei_LinearSystem.hpp>
//---------------------------------------------------------------------------
Solver_AztecOO::Solver_AztecOO()
  : tolerance_(1.e-6), maxIters_(500), useTranspose_(false), paramlist_(NULL),
    linProb(NULL),
    useML_(false),
#ifdef HAVE_FEI_ML
   ml_prec_(NULL), ml_defaults_set_(false),
   ml_aztec_options_(NULL), ml_aztec_params_(NULL),
#endif
  name_(),
  dbgprefix_("SlvAzoo: ")
{
  azoo_ = new AztecOO;
  paramlist_ = new Teuchos::ParameterList;

#ifdef HAVE_FEI_ML
  ml_aztec_options_ = new int[AZ_OPTIONS_SIZE];
  ml_aztec_params_ = new double[AZ_PARAMS_SIZE];
#endif
}

//---------------------------------------------------------------------------
Solver_AztecOO::~Solver_AztecOO()
{
  delete azoo_;
  delete paramlist_;
  delete linProb;
#ifdef HAVE_FEI_ML
  delete [] ml_aztec_options_;
  delete [] ml_aztec_params_;
  delete ml_prec_;
#endif

  AZ_manage_memory(0, AZ_CLEAR_ALL, 0, NULL, NULL);
}

//---------------------------------------------------------------------------
AztecOO& Solver_AztecOO::getAztecOO()
{
  return( *azoo_ );
}

//---------------------------------------------------------------------------
void Solver_AztecOO::setUseML(bool useml)
{
  useML_ = useml;
}

//---------------------------------------------------------------------------
Teuchos::ParameterList& Solver_AztecOO::get_ParameterList()
{
  if (paramlist_ == NULL) {
    paramlist_ = new Teuchos::ParameterList;
  }

  return( *paramlist_ );
}

//---------------------------------------------------------------------------
int Solver_AztecOO::solve(fei::LinearSystem* linearSystem,
			  fei::Matrix* preconditioningMatrix,
			  const fei::ParameterSet& parameterSet,
			  int& iterationsTaken,
			  int& status)
{
  std::string pcstring;
  parameterSet.getStringParamValue("AZ_precond", pcstring);
  if (pcstring == "ML_Op") {
    useML_ = true;
  }

  Teuchos::ParameterList& paramlist = get_ParameterList();

#ifdef HAVE_FEI_ML
  if (ml_aztec_options_ == NULL)
    ml_aztec_options_ = new int[AZ_OPTIONS_SIZE];
  if (ml_aztec_params_ == NULL)
    ml_aztec_params_ = new double[AZ_PARAMS_SIZE];

  if (!ml_defaults_set_ && useML_) {
    Teuchos::ParameterList mlparams;
    ML_Epetra::SetDefaults("SA", mlparams, ml_aztec_options_,ml_aztec_params_);
    mlparams.setParameters(paramlist);
    paramlist = mlparams;
    ml_defaults_set_ = true;
  }
#endif

  Trilinos_Helpers::copy_parameterset(parameterSet, paramlist);

  fei::SharedPtr<fei::Matrix> feiA = linearSystem->getMatrix();
  fei::SharedPtr<fei::Vector> feix = linearSystem->getSolutionVector();
  fei::SharedPtr<fei::Vector> feib = linearSystem->getRHS();

  Epetra_MultiVector*    x = NULL;
  Epetra_MultiVector*    b = NULL;
  Epetra_Operator* epetra_op = 0;
  Epetra_CrsMatrix* crsA = NULL;

  Trilinos_Helpers::get_Epetra_pointers(feiA, feix, feib,
                                        crsA, epetra_op, x, b);

  if (epetra_op == 0 || x == 0 || b == 0) {
    fei::console_out() << "Solver_AztecOO::solve Error, couldn't obtain Epetra objects"
     << " from fei container-objects."<<FEI_ENDL;
    return(-1);
  }

  //when we call azoo_->SetProblem, it will set some options. So we will
  //first take a copy of all options and params, then reset them after the
  //call to SetProblem. That way we preserve any options that have already
  //been set.

  std::vector<int> azoptions(AZ_OPTIONS_SIZE);
  std::vector<double> azparams(AZ_PARAMS_SIZE);

  const int* azoptionsptr = azoo_->GetAllAztecOptions();
  const double* azparamsptr = azoo_->GetAllAztecParams();

  int i;
  for(i=0; i<AZ_OPTIONS_SIZE; ++i) {
    azoptions[i] = azoptionsptr[i];
  }
  for(i=0; i<AZ_PARAMS_SIZE; ++i) {
    azparams[i] = azparamsptr[i];
  }

  Epetra_RowMatrix* precond = NULL;
  if (preconditioningMatrix != NULL) {
    fei::Matrix_Impl<Epetra_CrsMatrix>* snl_epetra_crs =
      dynamic_cast<fei::Matrix_Impl<Epetra_CrsMatrix>*>(preconditioningMatrix);
    fei::Matrix_Impl<Epetra_VbrMatrix>* snl_epetra_vbr =
      dynamic_cast<fei::Matrix_Impl<Epetra_VbrMatrix>*>(preconditioningMatrix);
    if (snl_epetra_crs != NULL) {
      precond = snl_epetra_crs->getMatrix().get();
    }
    else if (snl_epetra_vbr != NULL) {
      precond = snl_epetra_vbr->getMatrix().get();
    }
    else {
      fei::console_out() << "Solver_AztecOO::solve: ERROR getting epetra row matrix"
	       << " from preconditioningMatrix."<<FEI_ENDL;
      return(-1);
    }
  }

  if (precond != NULL) {
    Epetra_LinearProblem * newlinProb = new Epetra_LinearProblem(epetra_op,x,b);
    azoo_->SetProblem(*newlinProb);
    delete linProb;
    linProb = newlinProb;

    azoo_->SetAllAztecOptions(&(azoptions[0]));
    azoo_->SetAllAztecParams(&(azparams[0]));

    azoo_->SetUseAdaptiveDefaultsTrue();

    azoo_->SetPrecMatrix(precond);
  }

  bool needNewPreconditioner = false;

  if (feiA->changedSinceMark()) {
    feiA->markState();
    needNewPreconditioner = true;
  }

  if (needNewPreconditioner) {
    Epetra_LinearProblem * newlinProb = new Epetra_LinearProblem(epetra_op,x,b);
    azoo_->SetProblem(*newlinProb);
    delete linProb;
    linProb = newlinProb;

    azoo_->SetAllAztecOptions(&(azoptions[0]));
    azoo_->SetAllAztecParams(&(azparams[0]));

    azoo_->SetUseAdaptiveDefaultsTrue();

    if (useML_) {
#ifdef HAVE_FEI_ML
      setup_ml_operator(*azoo_, crsA);
#else
      fei::console_out() <<"Solver_AztecOO::solve ERROR, ML requested but HAVE_FEI_ML not defined."
	       << FEI_ENDL;
      return(-1);
#endif
    }
    else {
      azoo_->SetAztecOption(AZ_pre_calc, AZ_calc);
      azoo_->SetAztecOption(AZ_keep_info, 1);
    }
  }
  else {
    if (!useML_) {
      azoo_->SetAztecOption(AZ_pre_calc, AZ_reuse);
    }
  }

  epetra_op->SetUseTranspose(useTranspose_);

  azoo_->SetParameters(paramlist);

  maxIters_ = azoptionsptr[AZ_max_iter];
  tolerance_= azparamsptr[AZ_tol];

  status = azoo_->Iterate(maxIters_,
			 //2,//maxSolveAttempts
			 tolerance_);

  iterationsTaken = azoo_->NumIters();

  int olevel = 0;
  parameterSet.getIntParamValue("outputLevel", olevel);

  std::string param2;
  parameterSet.getStringParamValue("FEI_OUTPUT_LEVEL", param2);

  if (olevel >= 3 || param2 == "MATRIX_FILES" || param2 == "ALL") {
    std::string param1;
    parameterSet.getStringParamValue("debugOutput", param1);

    FEI_OSTRINGSTREAM osstr;
    if (!param1.empty()) {
      osstr << param1 << "/";
    }
    else osstr << "./";

    osstr << "x_AztecOO.vec";
    feix->writeToFile(osstr.str().c_str());
  }

  return(0);
}

//---------------------------------------------------------------------------
int Solver_AztecOO::solve(fei::LinearSystem* linearSystem,
			  fei::Matrix* preconditioningMatrix,
			  int numParams,
			  const char* const* solverParams,
			  int& iterationsTaken,
			  int& status)
{
  std::vector<std::string> stdstrings;
  fei::utils::char_ptrs_to_strings(numParams, solverParams, stdstrings);

  fei::ParameterSet paramset;
  fei::utils::parse_strings(stdstrings, " ", paramset);

  return( solve(linearSystem, preconditioningMatrix, paramset,
		iterationsTaken, status) );
}

//---------------------------------------------------------------------------
int Solver_AztecOO::setup_ml_operator(AztecOO& azoo, Epetra_CrsMatrix* A)
{
#ifdef HAVE_FEI_ML
  if (ml_aztec_options_ == NULL) {
    ml_aztec_options_ = new int[AZ_OPTIONS_SIZE];
  }
  if (ml_aztec_params_ == NULL) {
    ml_aztec_params_ = new double[AZ_PARAMS_SIZE];
  }

  if (!ml_defaults_set_) {
    Teuchos::ParameterList mlparams;
    ML_Epetra::SetDefaults("SA", mlparams,ml_aztec_options_,ml_aztec_params_);
    mlparams.setParameters(*paramlist_);
    *paramlist_ = mlparams;
    ml_defaults_set_ = true;
  }

  if (ml_prec_ != NULL) {
    delete ml_prec_; ml_prec_ = NULL;
  }

  ml_prec_ = new ML_Epetra::MultiLevelPreconditioner(*A, *paramlist_, true);
  azoo_->SetPrecOperator(ml_prec_);

#endif

  return(0);
}

#endif
//HAVE_FEI_AZTECOO

