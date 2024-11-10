// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "AztecOOParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_ValidatorXMLConverterDB.hpp"
#include "Teuchos_StandardValidatorXMLConverters.hpp"

namespace {

//
// Define the names of the different parameters.  Since the name of a
// parameter is used several times, it is a good idea to define a variable that
// stores the std::string name so that typing errors get caught at compile-time.
//

const std::string AztecSolver_name = "Aztec Solver";

const std::string AztecPreconditioner_name = "Aztec Preconditioner";

enum EAztecPreconditioner {
  AZTEC_PREC_NONE, AZTEC_PREC_ILU, AZTEC_PREC_ILUT, AZTEC_PREC_JACOBI,
  AZTEC_PREC_SYMMGS, AZTEC_PREC_POLY, AZTEC_PREC_LSPOLY
};

/** Needed for serialization KLN 23/09/2010 */
inline std::istream& operator>>(std::istream& is, EAztecPreconditioner& prec){
  int intval;
  is >> intval;
  prec = (EAztecPreconditioner)intval;
  return is;
}

const std::string Overlap_name = "Overlap";

const std::string GraphFill_name = "Graph Fill";

const std::string DropTolerance_name = "Drop Tolerance";

const std::string FillFactor_name = "Fill Factor";

const std::string Steps_name = "Steps";

const std::string PolynomialOrder_name = "Polynomial Order";

const std::string RCMReordering_name = "RCM Reordering";

const std::string Orthogonalization_name = "Orthogonalization";

const std::string SizeOfKrylovSubspace_name = "Size of Krylov Subspace";

const std::string ConvergenceTest_name = "Convergence Test";

const std::string IllConditioningThreshold_name = "Ill-Conditioning Threshold";

const std::string OutputFrequency_name = "Output Frequency";

Teuchos::RCP<Teuchos::ParameterList>  validAztecOOParams;

} // namespace

void setAztecOOParameters(
  Teuchos::ParameterList  *pl
  ,AztecOO                *solver
  )
{
  using Teuchos::getIntegralValue;
  using Teuchos::getParameter;
  TEUCHOS_TEST_FOR_EXCEPT(pl==NULL);
  TEUCHOS_TEST_FOR_EXCEPT(solver==NULL);
  // Validate the parameters and set their defaults!  This also sets the
  // validators needed to read in the parameters in an integral form.
  pl->validateParametersAndSetDefaults(*getValidAztecOOParameters());
  // Aztec Solver
  solver->SetAztecOption(
    AZ_solver
    ,getIntegralValue<int>(*pl,AztecSolver_name)
    );
  // Aztec Preconditioner
  switch(
    getIntegralValue<EAztecPreconditioner>(
      *pl,AztecPreconditioner_name
      )
    )
  {
    // This is the only place that EAztecPreconditioner is used.  Everywhere
    // else the code expects a string value.
    case AZTEC_PREC_NONE:
      solver->SetAztecOption(AZ_precond,AZ_none);
      pl->set(AztecPreconditioner_name, "none");
      break;
    case AZTEC_PREC_ILU:
      solver->SetAztecOption(AZ_precond,AZ_dom_decomp);
      solver->SetAztecOption(AZ_overlap,getParameter<int>(*pl,Overlap_name));
      solver->SetAztecOption(AZ_subdomain_solve,AZ_ilu);
      solver->SetAztecOption(AZ_graph_fill,getParameter<int>(*pl,GraphFill_name));
      pl->set(AztecPreconditioner_name, "ilu");
      break;
    case AZTEC_PREC_ILUT:
      solver->SetAztecOption(AZ_precond,AZ_dom_decomp);
      solver->SetAztecOption(AZ_overlap,getParameter<int>(*pl,Overlap_name));
      solver->SetAztecOption(AZ_subdomain_solve,AZ_ilut);
      solver->SetAztecParam(AZ_drop,getParameter<double>(*pl,DropTolerance_name));
      solver->SetAztecParam(AZ_ilut_fill,getParameter<double>(*pl,FillFactor_name));
      pl->set(AztecPreconditioner_name, "ilut");
      break;
    case AZTEC_PREC_JACOBI:
      solver->SetAztecOption(AZ_precond,AZ_Jacobi);
      solver->SetAztecOption(AZ_poly_ord,getParameter<int>(*pl,Steps_name));
      pl->set(AztecPreconditioner_name, "Jacobi");
      break;
    case AZTEC_PREC_SYMMGS:
      solver->SetAztecOption(AZ_precond,AZ_sym_GS);
      solver->SetAztecOption(AZ_poly_ord,getParameter<int>(*pl,Steps_name));
      pl->set(AztecPreconditioner_name, "Symmetric Gauss-Seidel");
      break;
    case AZTEC_PREC_POLY:
      solver->SetAztecOption(AZ_precond,AZ_Neumann);
      solver->SetAztecOption(AZ_poly_ord,getParameter<int>(*pl,PolynomialOrder_name));
      pl->set(AztecPreconditioner_name, "Polynomial");
      break;
    case AZTEC_PREC_LSPOLY:
      solver->SetAztecOption(AZ_precond,AZ_ls);
      solver->SetAztecOption(AZ_poly_ord,getParameter<int>(*pl,PolynomialOrder_name));
      pl->set(AztecPreconditioner_name, "Least-squares Polynomial");
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true); // Should never get here!
  }
  // RCM Reordering (in conjunction with domain decomp preconditioning)
  solver->SetAztecOption(
    AZ_reorder
    ,getIntegralValue<int>(*pl,RCMReordering_name)
    );
  // Gram-Schmidt orthogonalization procedure (GMRES only)
  solver->SetAztecOption(
    AZ_orthog
    ,getIntegralValue<int>(*pl,Orthogonalization_name)
    );
  // Size of the krylov subspace
  solver->SetAztecOption(AZ_kspace,getParameter<int>(*pl,SizeOfKrylovSubspace_name));
  // Convergence criteria to use in the linear solver
  solver->SetAztecOption(
    AZ_conv
    ,getIntegralValue<int>(*pl,ConvergenceTest_name)
    );
  // Set the ill-conditioning threshold for the upper hessenberg matrix
  solver->SetAztecParam(
    AZ_ill_cond_thresh, getParameter<double>(*pl,IllConditioningThreshold_name)
    );
  // Frequency of linear solve residual output
  solver->SetAztecOption(
    AZ_output, getParameter<int>(*pl,OutputFrequency_name)
    );
#ifdef TEUCHOS_DEBUG
  // Check to make sure that I did not use the PL incorrectly!
  pl->validateParameters(*getValidAztecOOParameters());
#endif // TEUCHOS_DEBUG
}

Teuchos::RCP<const Teuchos::ParameterList>
getValidAztecOOParameters()
{
  //
  // This function defines the valid parameter list complete with validators
  // and default values.  The default values never need to be repeated because
  // if the use of the function validateParametersAndSetDefaults(...) used
  // above in setAztecOOParameters(...).  Also, the validators do not need to
  // be kept track of since they will be set in the input list also.
  //
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::setIntParameter;
  using Teuchos::setDoubleParameter;
  using Teuchos::ParameterList;
  /**
   * Must add this otherwise one of the stratimikos tests fail.
   * By adding this converter we can write the ParameterList out to XML
   * successfully.
   * KLN 23/09/2010
   */
  Teuchos::ValidatorXMLConverterDB::addConverter(
    Teuchos::DummyObjectGetter<
      Teuchos::StringToIntegralParameterEntryValidator<EAztecPreconditioner> 
    >::getDummyObject(),
    Teuchos::DummyObjectGetter<Teuchos::StringToIntegralValidatorXMLConverter<
      EAztecPreconditioner> >::getDummyObject());
  //
  RCP<ParameterList> pl = validAztecOOParams;
  if(pl.get()) return pl;
  pl = validAztecOOParams = rcp(new ParameterList());
  //
  setStringToIntegralParameter<int>(
    AztecSolver_name, "GMRES",
    "Type of linear solver algorithm to use.",
    tuple<std::string>("CG","GMRES","CGS","TFQMR","BiCGStab","LU","GMRESR","FixedPoint"),
    tuple<int>(AZ_cg,AZ_gmres,AZ_cgs,AZ_tfqmr,AZ_bicgstab,AZ_lu,AZ_GMRESR,AZ_fixed_pt),
    &*pl
    );
  setStringToIntegralParameter<EAztecPreconditioner>(
    AztecPreconditioner_name, "ilu",
    "Type of internal preconditioner to use.\n"
    "Note! this preconditioner will only be used if the input operator\n"
    "supports the Epetra_RowMatrix interface and the client does not pass\n"
    "in an external preconditioner!",
    tuple<std::string>(
      "none","ilu","ilut","Jacobi",
      "Symmetric Gauss-Seidel","Polynomial","Least-squares Polynomial"
      ),
    tuple<EAztecPreconditioner>(
      AZTEC_PREC_NONE,AZTEC_PREC_ILU,AZTEC_PREC_ILUT,AZTEC_PREC_JACOBI,
      AZTEC_PREC_SYMMGS,AZTEC_PREC_POLY,AZTEC_PREC_LSPOLY
      ),
    &*pl
    );
  setIntParameter(
    Overlap_name, 0,
    "The amount of overlap used for the internal \"ilu\" and \"ilut\" preconditioners.",
    &*pl
    );
  setIntParameter(
    GraphFill_name, 0,
    "The amount of fill allowed for the internal \"ilu\" preconditioner.",
    &*pl
    );
  setDoubleParameter(
    DropTolerance_name, 0.0,
    "The tolerance below which an entry from the factors of an internal \"ilut\"\n"
    "preconditioner will be dropped.",
    &*pl
    );
  setDoubleParameter(
    FillFactor_name, 1.0,
    "The amount of fill allowed for an internal \"ilut\" preconditioner.",
    &*pl
    );
  setIntParameter(
    Steps_name, 3,
    "Number of steps taken for the \"Jacobi\" or the \"Symmetric Gauss-Seidel\"\n"
    "internal preconditioners for each preconditioner application.",
    &*pl
    );
  setIntParameter(
    PolynomialOrder_name, 3,
    "The order for of the polynomials used for the \"Polynomial\" and\n"
    "\"Least-squares Polynomial\" internal preconditioners.",
    &*pl
    );
  setStringToIntegralParameter<int>(
    RCMReordering_name, "Disabled",
    "Determines if RCM reordering is used with the internal\n"
    "\"ilu\" or \"ilut\" preconditioners.",
    tuple<std::string>("Enabled","Disabled"),
    tuple<int>(1,0),
    &*pl
    );
  setStringToIntegralParameter<int>(
    Orthogonalization_name, "Classical",
    "The type of orthogonalization to use with the \"GMRES\" solver.",
    tuple<std::string>("Classical","Modified"),
    tuple<int>(AZ_classic,AZ_modified),
    &*pl
    );
  setIntParameter(
    SizeOfKrylovSubspace_name, 300,
    "The maximum size of the Krylov subspace used with \"GMRES\" before\n"
    "a restart is performed.",
    &*pl
    );
  setStringToIntegralParameter<int>(
    ConvergenceTest_name, "r0", // Same as "rhs" when x=0
    "The convergence test to use for terminating the iterative solver.",
    tuple<std::string>("r0","rhs","Anorm","no scaling","sol"),
    tuple<int>(AZ_r0,AZ_rhs,AZ_Anorm,AZ_noscaled,AZ_sol),
    &*pl
    );
  setDoubleParameter(
    IllConditioningThreshold_name, 1e+11,
    "The threshold tolerance above which a system is considered\n"
    "ill conditioned.",
    &*pl
    );
  setIntParameter(
    OutputFrequency_name, 0, // By default, no output from Aztec!
    "The number of iterations between each output of the solver's progress.",
    &*pl
    );
  //
  return pl;
}
