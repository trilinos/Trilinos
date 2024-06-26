// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
// This driver tests the translation of Aztec parameters to 
// Teuchos::ParameterList based Belos compatable parameters. 
//
// NOTE: No preconditioner or subdomain solve translateion is supported, aand
// will generate a Belos::XLateStatus::WARN as part of the return pair;
// 
//
#include "BelosConfigDefs.hpp"
#include "BelosUtils.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "az_aztec_defs.h"

using namespace Teuchos;

int main(int argc, char *argv[]) {
  
  //* indicates it is a option that is translated by Belos_translateFromAztecParams
  //x indicates that this option will produce a warning return type. 
//AZ_solver              0   *
//AZ_scaling             1
//AZ_precond             2 x
//AZ_conv                3   *
//AZ_output              4   *
//AZ_pre_calc            5
//AZ_max_iter            6   *
//AZ_poly_ord            7
//AZ_overlap             8
//AZ_type_overlap        9
//AZ_kspace              10   * if AZ_solver == AZ_gmres
//AZ_orthog              11   * ICGS and IMGS only
//AZ_aux_vec             12
//AZ_reorder             13
//AZ_keep_info           14
//AZ_recursion_level     15
//AZ_print_freq          16
//AZ_graph_fill          17
//AZ_subdomain_solve     18 x
//AZ_init_guess          19
//AZ_keep_kvecs          20
//AZ_apply_kvecs         21
//AZ_orth_kvecs          22
//AZ_ignore_scaling      23
//AZ_check_update_size   24
//AZ_extreme             25
//AZ_diagnostics         26

// AZ_tol                 0   *
// AZ_drop                1
// AZ_ilut_fill           2
// AZ_omega               3
// AZ_rthresh             4
// AZ_athresh             5
// AZ_update_reduction    6
// AZ_temp                7
// AZ_ill_cond_thresh     8
// AZ_weights             9

 std::pair<std::string,int> xlate_err;

  int * az_opt = new int [AZ_FIRST_USER_OPTION];
  memset(az_opt,0,sizeof(int)*(AZ_FIRST_USER_OPTION));
  std::vector<double> Vaz_param;
  Vaz_param.resize(AZ_FIRST_USER_PARAM,0.0);

  Vaz_param[AZ_tol] = 1e-9;
  // treatment of specified solvers. 
// AZ_cg               0  *
// AZ_gmres            1  *
// AZ_cgs              2  x Error
// AZ_tfqmr            3  *
// AZ_bicgstab         4  *
// AZ_slu              5  x Error
// AZ_symmlq           6  x Error
// AZ_GMRESR           7  x Error
// AZ_fixed_pt         8  x Error
// AZ_analyze          9  x Error
// AZ_lu              10  *
// AZ_cg_condnum      11  *
// AZ_gmres_condnum   12  x Error

  az_opt[AZ_solver] = AZ_gmres;
  az_opt[AZ_conv] = AZ_r0; // or AZ_rhs or AZ_Anorm or AZ_noscaled
  
  Teuchos::ParameterList tpl;
  const double * az_par_val = (const double * ) &(Vaz_param[0]);
  xlate_err = Belos::translateFromAztecParams(tpl,az_opt,az_par_val);
  
  if(xlate_err.second != Belos::TRANSLATE_FROM_AZTEC_OK || xlate_err.first.size()!=0 ) {
    // this one should be error and warning free
    std::cout << " translateFromAztecParams::  failure, string is:"<<std::endl;
    std::cout << xlate_err.first<<std::endl;
    std::cout << " Error num "<< xlate_err.second<<std::endl;
    std::cout << " enum   Opt Param "<<std::endl;
    for(int i=0;i<AZ_FIRST_USER_OPTION;++i) {
      std::cout << i<<" "<<az_opt[i];
      if(i<AZ_FIRST_USER_PARAM) std::cout <<" "<<az_par_val[i];
      std::cout<<std::endl;
    }
   
    tpl.print();
    return EXIT_FAILURE;
  }
  tpl.print();
  std::cout<<" Pass "<<std::endl;

  // now add some stuff that should give warnings.
  az_opt[AZ_precond] = AZ_ls; // any value other than AZ_none == 0 generates a warning. 
  az_opt[AZ_conv] = AZ_rhs; // should be valid. 
  
  xlate_err = Belos::translateFromAztecParams(tpl,az_opt,az_par_val);
  
  if(xlate_err.second != Belos::TRANSLATE_FROM_AZTEC_WARN) {
    // this one should be error free but generate a  warning
    std::cout << " translateFromAztecParams::  failure, string is:"<<std::endl;
    std::cout << xlate_err.first<<std::endl;
    std::cout << " Error num "<< xlate_err.second<<std::endl;
    tpl.print();
    return EXIT_FAILURE;
  }

  tpl.print();
  std::cout<<" Pass "<<std::endl;
  
  az_opt[AZ_precond] = AZ_none;
  az_opt[AZ_subdomain_solve] = AZ_icc;
  az_opt[AZ_conv] = AZ_Anorm;

  xlate_err = Belos::translateFromAztecParams(tpl,az_opt,az_par_val);
  
  if(xlate_err.second != Belos::TRANSLATE_FROM_AZTEC_WARN) {
    // this one should be error free but generate a  warning
    std::cout << " translateFromAztecParams::  failure, string is:"<<std::endl;
    std::cout << xlate_err.first<<std::endl;
    std::cout << " Error num "<< xlate_err.second<<std::endl;
    tpl.print();
    return EXIT_FAILURE;
  }

  tpl.print();
  std::cout<<" Pass "<<std::endl;

  // now errors
  az_opt[AZ_orthog]=AZ_double_classic;
  az_opt[AZ_conv] = AZ_noscaled;
  az_opt[AZ_subdomain_solve] = AZ_lu;
  xlate_err = Belos::translateFromAztecParams(tpl,az_opt,az_par_val);
  
  if(! (xlate_err.second | Belos::TRANSLATE_FROM_AZTEC_ERROR && xlate_err.second|Belos::TRANSLATE_FROM_AZTEC_WARN)) {
    // This should generate an error and a warning. 
    // error from az_double_classic
    std::cout << " translateFromAztecParams::  failure, string is:"<<std::endl;
    std::cout << xlate_err.first<<std::endl;
    std::cout << " Error num "<< xlate_err.second<<std::endl;
    tpl.print();
    return EXIT_FAILURE;
  }

  tpl.print();
  std::cout<<" Pass All"<<std::endl;
  return EXIT_SUCCESS;
}
  
  
  
