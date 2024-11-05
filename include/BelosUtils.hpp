// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_UTIL_HPP
#define BELOS_UTIL_HPP
#include "BelosConfigDefs.hpp"

#ifdef HAVE_BELOS_AZTECOO
#include "az_aztec_defs.h" 
#include "Teuchos_ParameterList.hpp"
#include "BelosTypes.hpp"

namespace Belos{
 enum ETranslateFromAztecStatus {
    TRANSLATE_FROM_AZTEC_OK    =0x00,
    TRANSLATE_FROM_AZTEC_WARN  =0x01,
    TRANSLATE_FROM_AZTEC_ERROR =0x02};

std::pair< std::string, int >  
translateFromAztecParams( Teuchos::ParameterList &tpl ,
			  const int * aztec_options, 
			  const double * aztec_params
			  ) {
  using namespace std;
  int econd = TRANSLATE_FROM_AZTEC_OK;
  ostringstream error;
  if(aztec_options == NULL || aztec_params == NULL ) {
    return std::pair<std::string,int>(string("Belos_Translate_from_Aztec_Params:: Aztec Options or Parameters were null."),econd);
  }

  switch (aztec_options[AZ_solver]){
  case AZ_gmres:

    tpl.set("Solver Name","Pseudoblock GMRES");
    break;
  case AZ_cg:
  case AZ_cg_condnum:
    tpl.set("Solver Name","Pseudoblock CG");
    break; 
  case AZ_bicgstab:
    tpl.set("Solver Name","BICGSTAB");
    break;
  case AZ_tfqmr:
    tpl.set("Solver Name","TFQMR");
    break;
  case AZ_fixed_pt:
    tpl.set("Solver Name","FIXED POINT");
    break;
  case AZ_lu:
    error<<" Translate_Params_Aztec_to_Belos:: uncaught solver name  AZ_lu "<<std::endl;
    econd |= TRANSLATE_FROM_AZTEC_ERROR;
    break;
  case AZ_cgs:
    error<<" Translate_Params_Aztec_to_Belos:: uncaught solver name  AZ_cgs"<<std::endl;
    econd |= TRANSLATE_FROM_AZTEC_ERROR;
    break;
  case AZ_slu:
    error<<" Translate_Params_Aztec_to_Belos:: uncaught solver name  AZ_slu"<<std::endl;
    econd |= TRANSLATE_FROM_AZTEC_ERROR;
    break;
  case AZ_symmlq:
    error<<" Translate_Params_Aztec_to_Belos:: uncaught solver name  AZ_symmlq"<<std::endl;
    econd |= TRANSLATE_FROM_AZTEC_ERROR;
    break;
  case AZ_GMRESR:
    error<<" Translate_Params_Aztec_to_Belos:: uncaught solver name  AZ_GMRESR"<<std::endl;
    econd |= TRANSLATE_FROM_AZTEC_ERROR;
    break;

  case AZ_analyze:
    error<<" Translate_Params_Aztec_to_Belos:: uncaught solver name  AZ_analyze "<<std::endl;
    econd |= TRANSLATE_FROM_AZTEC_ERROR;
    break;
  case AZ_gmres_condnum:
    error<<" Translate_Params_Aztec_to_Belos:: uncaught solver name  AZ_gmres_condnum."<<std::endl;
    econd |= TRANSLATE_FROM_AZTEC_ERROR;
    break;

  default:
    error << "Translate_Params_Aztec_to_Belos:: unknown solver enum "<<aztec_options[AZ_solver]<<std::endl;
    econd |= TRANSLATE_FROM_AZTEC_ERROR;
  }
  // sierra
  //PRECONDITIONING METHOD=DD-ICC
  //PRECONDITIONING METHOD=DD-ILU
  //PRECONDITIONING METHOD=DD-ILUT

  switch (aztec_options[AZ_precond]) {
  case AZ_none:
    break; // only valid option.
  case AZ_sym_GS: 
  case AZ_ls:
  case AZ_Jacobi:
  case AZ_Neumann:
  case AZ_dom_decomp:
  default:
    error<<" Belos does not have built in preconditioners, Az_precond ignored."<<std::endl;
    econd |= TRANSLATE_FROM_AZTEC_WARN;
  };

  switch(aztec_options[AZ_subdomain_solve]) {
  case AZ_none:
    break; // only valid option
  case AZ_lu:
  case AZ_ilut:
  case AZ_ilu:
  case AZ_rilu:
  case AZ_bilu:
  case AZ_icc:
  default:
      error<<" Belos does not have built in subdomain solvers, Az_subdomain_solve ignored."<<std::endl;
    econd |= TRANSLATE_FROM_AZTEC_WARN;
  }

  // All sierra options are here.
  switch(aztec_options[AZ_conv]) {
  case AZ_r0:
    tpl.set("Implicit Residual Scaling","Norm of Initial Residual");
    break;
  case AZ_rhs:
    tpl.set("Implicit Residual Scaling","Norm of RHS");
    break;
  case AZ_Anorm:
    tpl.set("Implicit Residual Scaling","Norm of Preconditioned Initial Residual");
    break;
  case AZ_noscaled:
    tpl.set("Implicit Residual Scaling","None");
    break;
  case AZ_sol:   
  case AZ_expected_values:
  default: 
    error << "Belos_Translate_from_Aztec_Params: AZ_conv of AZ_sol or AZ_expected_values are not valid for belos. "<<std::endl;
    econd |= TRANSLATE_FROM_AZTEC_ERROR;
  }

  // Make Belos produce output like AztecOO's.
  tpl.set("Output Style", static_cast<int> (Belos::Brief));
  // Always print Belos' errors.  You can add the 
  // enum values together to get all of their effects.
  Belos::MsgType belosPrintOptions = static_cast<Belos::MsgType> (static_cast<int> (Belos::Errors) + static_cast<int> (Belos::Warnings));

  switch(aztec_options[AZ_output]) {
    //  "Errors",
    //     "Warnings",
    //     "IterationDetails",
    //     "OrthoDetails",
    //     "FinalSummary",
    //     "TimingDetails",
    //     "StatusTestDetails",
    //     "Debug"
    
  case AZ_none:
    tpl.set("Output Frequency", -1);
    break;
  case AZ_all:
    tpl.set("Output Frequency", 1);
    belosPrintOptions = static_cast<Belos::MsgType> (static_cast<int> (belosPrintOptions) + static_cast<int> (Belos::StatusTestDetails));
    belosPrintOptions = static_cast<Belos::MsgType> (static_cast<int> (belosPrintOptions) + static_cast<int> (Belos::FinalSummary));
    break;
  case AZ_warnings: // only print warnings
    tpl.set("Output Frequency", -1);
    break;
  case AZ_last:	// only print the final result
    tpl.set("Output Frequency", -1);
    belosPrintOptions = static_cast<Belos::MsgType> (static_cast<int> (belosPrintOptions) + static_cast<int> (Belos::FinalSummary));
    break;
  default: // some integer; print every that many iterations
    // This is an int and not a enum, so I don't need to cast it.
    const int freq = aztec_options[AZ_output];
    tpl.set("Output Frequency", freq);
    belosPrintOptions = static_cast<Belos::MsgType> (static_cast<int> (belosPrintOptions) + static_cast<int> (Belos::StatusTestDetails));
    belosPrintOptions = static_cast<Belos::MsgType> (static_cast<int> (belosPrintOptions) + static_cast<int> (Belos::FinalSummary));
    break;
  }
  tpl.set("Verbosity", static_cast<int> (belosPrintOptions));

  // Only set the "Orthogonalization" parameter if using GMRES.  CG
  // doesn't accept that parameter.
  // sierra uses only ORTHOG METHOD=CLASSICAL
  if (aztec_options[AZ_solver] == AZ_gmres) {
    switch(aztec_options[AZ_orthog]) {
    case AZ_classic:
      tpl.set("Orthogonalization","ICGS");
      break;
    case AZ_modified:
      tpl.set("Orthogonalization","IMGS");
      break;
    default:
      error<<"Option AZ_orthog for GMRES not recognized "<<aztec_options[AZ_orthog]<<endl;
      econd |= TRANSLATE_FROM_AZTEC_ERROR;
      // no way to set DGKS
    }
  }

  if(aztec_options[AZ_max_iter]!=0)   
    tpl.set("Maximum Iterations",aztec_options[AZ_max_iter]);

 // Only set the "Num Blocks" (restart length) parameter if using
  // GMRES.  CG doesn't accept that parameter.
  if (aztec_options[AZ_solver] == AZ_gmres && 
      aztec_options[AZ_kspace] !=0) {
    tpl.set("Num Blocks",aztec_options[AZ_kspace]);
  }

  // aztec_params tested, only AZ_tol should be set. 
  tpl.set("Convergence Tolerance",aztec_params[AZ_tol]);
  for(int i=AZ_drop ; i<= AZ_weights ; ++i) {
    if(aztec_params[i]!=0 ){
      error << " Aztec_Params at "<<i<<" non zero and will be ignored"<<std::endl;
      econd |= TRANSLATE_FROM_AZTEC_WARN;
    }
  }

 

  return std::pair<std::string,int>(error.str(),econd);
}
}
#endif
#endif
