/*@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER
*/

#include "AztecOOParameterList.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StringToIntMap.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

namespace {

//
// Define the names of the different parameters, their default values
// and a validator if appropriate!
//
// All of this setup makes these parameters all validated and all
// maintainable!
//

const std::string AztecSolver_name = "Aztec Solver";
const Teuchos::RefCountPtr<Teuchos::StringToIntegralParameterEntryValidator<int> >
aztecSolverValidator = Teuchos::rcp(
  new Teuchos::StringToIntegralParameterEntryValidator<int>(
    Teuchos::tuple<std::string>("CG","GMRES","CGS","TFQMR","BiCGStab","LU")
    ,Teuchos::tuple<int>(AZ_cg,AZ_gmres,AZ_cgs,AZ_tfqmr,AZ_bicgstab,AZ_lu)
    ,AztecSolver_name
    )
  );
const std::string AztecSolver_default = "GMRES";

const std::string AztecPreconditioner_name = "Aztec Preconditioner";enum EAztecPreconditioner { AZTEC_PREC_NONE, AZTEC_PREC_ILU, AZTEC_PREC_ILUT, AZTEC_PREC_JACOBI, AZTEC_PREC_SYMMGS, AZTEC_PREC_POLY, AZTEC_PREC_LSPOLY };
const Teuchos::RefCountPtr<Teuchos::StringToIntegralParameterEntryValidator<EAztecPreconditioner> >
aztecPrecValidator = rcp(
  new Teuchos::StringToIntegralParameterEntryValidator<EAztecPreconditioner>(
    Teuchos::tuple<std::string>(
      "none","ilu","ilut","Jacobi","Symmetric Gauss-Seidel","Polynomial","Least-squares Polynomial"
      )
    ,Teuchos::tuple<EAztecPreconditioner>(
      AZTEC_PREC_NONE,AZTEC_PREC_ILU,AZTEC_PREC_ILUT,AZTEC_PREC_JACOBI,AZTEC_PREC_SYMMGS,AZTEC_PREC_POLY,AZTEC_PREC_LSPOLY
      )
    ,AztecPreconditioner_name
    )
  );
const std::string AztecPreconditioner_default = "ilu";
  
const std::string Overlap_name = "Overlap";
const int Overlap_default = 0;

const std::string GraphFill_name = "Graph Fill";
const int GraphFill_default = 0;

const std::string DropTolerance_name = "Drop Tolerance";
const double DropTolerance_default = 0.0;

const std::string FillFactor_name = "Fill Factor";
const double FillFactor_default = 1.0;

const std::string Steps_name = "Steps";
const int Steps_default = 3;

const std::string PolynomialOrder_name = "Polynomial Order";
const int PolynomialOrder_default = 3;

const std::string RCMReordering_name = "RCM Reordering";
const Teuchos::RefCountPtr<Teuchos::StringToIntegralParameterEntryValidator<int> >
rcmReorderingValidator = Teuchos::rcp(
  new Teuchos::StringToIntegralParameterEntryValidator<int>(
    Teuchos::tuple<std::string>("Enabled","Disabled")
    ,Teuchos::tuple<int>(1,0)
    ,RCMReordering_name
    )
  );
const std::string RCMReordering_default = "Disabled";

const std::string Orthogonalization_name = "Orthogonalization";
const Teuchos::RefCountPtr<Teuchos::StringToIntegralParameterEntryValidator<int> >
orthogValidator = Teuchos::rcp(
  new Teuchos::StringToIntegralParameterEntryValidator<int>(
    Teuchos::tuple<std::string>("Classical","Modified")
    ,Teuchos::tuple<int>(AZ_classic,AZ_modified)
    ,Orthogonalization_name
    )
  );
const std::string Orthogonalization_default = "Classical";

const std::string SizeOfKrylovSubspace_name = "Size of Krylov Subspace";
const int SizeOfKrylovSubspace_default = 300;

const std::string ConvergenceTest_name = "Convergence Test";
const Teuchos::RefCountPtr<Teuchos::StringToIntegralParameterEntryValidator<int> >
convTestValidator = Teuchos::rcp(
  new Teuchos::StringToIntegralParameterEntryValidator<int>(
    Teuchos::tuple<std::string>("r0","rhs","Anorm","no scaling","sol")
    ,Teuchos::tuple<int>(AZ_r0,AZ_rhs,AZ_Anorm,AZ_noscaled,AZ_sol)
    ,ConvergenceTest_name
    )
  );
const std::string ConvergenceTest_default = "r0"; // Same as "rhs" when x=0

const std::string IllConditioningThreshold_name = "Ill-Conditioning Threshold";
const double IllConditioningThreshold_default = 1e+11;

const std::string OutputFrequency_name = "Output Frequency";
const int OutputFrequency_default = 0; // Aztec only prints to stdout which is not acceptable!

Teuchos::RefCountPtr<Teuchos::ParameterList>  validAztecOOParams;

} // namespace

void setAztecOOParameters(
  Teuchos::ParameterList  *pl
  ,AztecOO                *solver
  )
{
  TEST_FOR_EXCEPT(pl==NULL);
  TEST_FOR_EXCEPT(solver==NULL);
  // Aztec Solver
  solver->SetAztecOption(
    AZ_solver
    ,aztecSolverValidator->getIntegralValue(*pl,AztecSolver_name,AztecSolver_default)
    );
  // Aztec Preconditioner
  switch(aztecPrecValidator->getIntegralValue(*pl,AztecPreconditioner_name,AztecPreconditioner_default))
  {
    case AZTEC_PREC_NONE:
      solver->SetAztecOption(AZ_precond,AZ_none);
      break;
    case AZTEC_PREC_ILU:
      solver->SetAztecOption(AZ_precond,AZ_dom_decomp);
      solver->SetAztecOption(AZ_overlap,pl->get(Overlap_name,Overlap_default));
      solver->SetAztecOption(AZ_subdomain_solve,AZ_ilu);
      solver->SetAztecOption(AZ_graph_fill,pl->get(GraphFill_name,GraphFill_default));
      break;
    case AZTEC_PREC_ILUT:
      solver->SetAztecOption(AZ_precond,AZ_dom_decomp);
      solver->SetAztecOption(AZ_overlap,pl->get(Overlap_name,Overlap_default));
      solver->SetAztecOption(AZ_subdomain_solve,AZ_ilut);
      solver->SetAztecParam(AZ_drop,pl->get(DropTolerance_name,DropTolerance_default));
      solver->SetAztecParam(AZ_ilut_fill,pl->get(FillFactor_name,FillFactor_default));
      break;
    case AZTEC_PREC_JACOBI:
      solver->SetAztecOption(AZ_precond,AZ_Jacobi);
      solver->SetAztecOption(AZ_poly_ord,pl->get(Steps_name,Steps_default));
      break;
    case AZTEC_PREC_SYMMGS:
      solver->SetAztecOption(AZ_precond,AZ_sym_GS);
      solver->SetAztecOption(AZ_poly_ord,pl->get(Steps_name,Steps_default));
      break;
    case AZTEC_PREC_POLY:
      solver->SetAztecOption(AZ_precond,AZ_Neumann);
      solver->SetAztecOption(AZ_poly_ord,pl->get(PolynomialOrder_name,PolynomialOrder_default));
      break;
    case AZTEC_PREC_LSPOLY:
      solver->SetAztecOption(AZ_precond,AZ_ls);
      solver->SetAztecOption(AZ_poly_ord,pl->get(PolynomialOrder_name,PolynomialOrder_default));
      break;
    default:
      TEST_FOR_EXCEPT(true); // Should never get here!
  }
  // RCM Reordering (in conjunction with domain decomp preconditioning)
  solver->SetAztecOption(
    AZ_reorder
    ,rcmReorderingValidator->getIntegralValue(*pl,RCMReordering_name,RCMReordering_default)
    );
  // Gram-Schmidt orthogonalization procedure (GMRES only)
  solver->SetAztecOption(
    AZ_orthog
    ,orthogValidator->getIntegralValue(*pl,Orthogonalization_name,Orthogonalization_default)
    );
  // Size of the krylov subspace
  solver->SetAztecOption(AZ_kspace,pl->get(SizeOfKrylovSubspace_name,SizeOfKrylovSubspace_default));
  // Convergence criteria to use in the linear solver
  solver->SetAztecOption(
    AZ_conv
    ,convTestValidator->getIntegralValue(*pl,ConvergenceTest_name,ConvergenceTest_default)
    );
  // Set the ill-conditioning threshold for the upper hessenberg matrix
  solver->SetAztecParam(
    AZ_ill_cond_thresh
    ,pl->get(IllConditioningThreshold_name,IllConditioningThreshold_default)
    );
  // Frequency of linear solve residual output
  solver->SetAztecOption(
    AZ_output
    ,pl->get(OutputFrequency_name,OutputFrequency_default)
    );
#ifdef TEUCHOS_DEBUG
  // Check to make sure that I did not use the PL incorrectly!
  pl->validateParameters(*getValidAztecOOParameters());
#endif // TEUCHOS_DEBUG
}

Teuchos::RefCountPtr<const Teuchos::ParameterList>
getValidAztecOOParameters()
{
  //
  using Teuchos::RefCountPtr;
  using Teuchos::rcp;
  using Teuchos::tuple;
  using Teuchos::StringToIntegralParameterEntryValidator;
  using Teuchos::ParameterList;
  //
  RefCountPtr<ParameterList> pl = validAztecOOParams;
  if(pl.get()) return pl;
  pl = validAztecOOParams = rcp(new ParameterList());
  //
  pl->set(
    AztecSolver_name,AztecSolver_default
    ,"Type of linear solver algorithm to use."
    ,aztecSolverValidator
    );
  pl->set(
    AztecPreconditioner_name,AztecPreconditioner_default
    ,"Type of internal preconditioner to use.\n"
    "Note! this preconditioner will only be used if the input operator\n"
    "supports the Epetra_RowMatrix interface and the client does not pass\n"
    "in an external preconditioner!"
    ,aztecPrecValidator
    );
  pl->set(
    Overlap_name,Overlap_default
    ,"The amount of overlap used for the internal \"ilu\" and \"ilut\" preconditioners."
    );
  pl->set(
    GraphFill_name,GraphFill_default
    ,"The amount of fill allowed for the internal \"ilu\" preconditioner."
    );
  pl->set(
    DropTolerance_name,DropTolerance_default
    ,"The tolerance below which an entry from the factors of an internal \"ilut\"\n"
    "preconditioner will be dropped."
    );
  pl->set(
    FillFactor_name,FillFactor_default
    ,"The amount of fill allowed for an internal \"ilut\" preconditioner."
    );
  pl->set(
    Steps_name,Steps_default
    ,"Number of steps taken for the \"Jacobi\" or the \"Symmetric Gauss-Seidel\"\n"
    "internal preconditioners for each preconditioner application."
    );
  pl->set(
    PolynomialOrder_name,PolynomialOrder_default
    ,"The order for of the polynomials used for the \"Polynomial\" and\n"
    "\"Least-squares Polynomial\" internal preconditioners."
    );
  pl->set(
    RCMReordering_name,"Disabled"
    ,"Determines if RCM reordering is used with the internal\n"
    "\"ilu\" or \"ilut\" preconditioners."
    ,rcmReorderingValidator
    );
  pl->set(
    Orthogonalization_name,Orthogonalization_default
    ,"The type of orthogonalization to use with the \"GMRES\" solver."
    ,orthogValidator
    );
  pl->set(
    SizeOfKrylovSubspace_name,SizeOfKrylovSubspace_default
    ,"The maximum size of the Krylov subspace used with \"GMRES\" before\n"
    "a restart is performed."
    );
  pl->set(
    ConvergenceTest_name,ConvergenceTest_default
    ,"The convergence test to use for terminating the iterative solver."
    ,convTestValidator
    );
  pl->set(
    IllConditioningThreshold_name,IllConditioningThreshold_default
    ,"The threshold tolerance above which a system is considered\n"
    "ill conditioned."
    );
  pl->set(
    OutputFrequency_name,OutputFrequency_default
    ,"The number of iterations between each output of the solver's progress."
    );
  //
  return pl;
}
