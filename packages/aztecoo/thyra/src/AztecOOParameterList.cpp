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

namespace {

//
// (Single defintion) parameter names
//

const std::string AztecSolver_name = "Aztec Solver";
const std::string AztecPreconditioner_name = "Aztec Preconditioner";
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

//
// Parameter value lists and (single definition) default values
//

// AztecSolver
const int numAztecSolverValues = 6;
enum EAztecSolver
{ AZTEC_SOLVER_CG, AZTEC_SOLVER_GMRES, AZTEC_SOLVER_CGS, AZTEC_SOLVER_TFQMR, AZTEC_SOLVER_BICGSTAB, AZTEC_SOLVER_LU };
const char* AztecSolverValues[numAztecSolverValues] =
{ "CG", "GMRES", "CGS", "TFQMR", "BiCGStab", "LU" };
const Teuchos::StringToIntMap
AztecSolverMap(AztecSolver_name,numAztecSolverValues,AztecSolverValues);
const EAztecSolver AztecSolver_default = AZTEC_SOLVER_GMRES;
// AztecPreconditioner
const int numAztecPreconditionerValues = 7;
enum EAztecPreconditioner
{ AZTEC_PREC_NONE, AZTEC_PREC_ILU, AZTEC_PREC_ILUT, AZTEC_PREC_JACOBI, AZTEC_PREC_SYMMGS, AZTEC_PREC_POLY, AZTEC_PREC_LSPOLY };
const char* AztecPreconditionerValues[numAztecPreconditionerValues] =
{ "none", "ilu", "ilut", "Jacobi", "Symmetric Gauss-Seidel", "Polynomial", "Least-squares Polynomial" };
const Teuchos::StringToIntMap
AztecPreconditionerMap(AztecPreconditioner_name,numAztecPreconditionerValues,AztecPreconditionerValues);
const EAztecPreconditioner AztecPreconditioner_default = AZTEC_PREC_ILU;
// Overlap
const int Overlap_default = 0;
// GraphFill
const int GraphFill_default = 0;
// DropTolerance
const double DropTolerance_default = 0.0;
// FillFactor
const double FillFactor_default = 1.0;
// Steps
const int Steps_default = 3;
// PolynomialOrder
const int PolynomialOrder_default = 3;
// RCMReordering
const int numRCMReorderingValues = 2;
enum ERCMReordering { RCM_REORDERING_ENABLED, RCM_REORDERING_DISABLED }; 
const char* RCMReorderingValues[numRCMReorderingValues] = { "Enabled", "Disabled" };
const Teuchos::StringToIntMap RCMReorderingMap(RCMReordering_name,numRCMReorderingValues,RCMReorderingValues);
const ERCMReordering RCMReordering_default = RCM_REORDERING_DISABLED;
// Orthogonalization
const int numOrthogonalizationValues = 2;
enum EOrthogonalization { ORTHOG_CLASSICAL, ORTHOG_MODIFIED }; 
const char* OrthogonalizationValues[numOrthogonalizationValues] = { "Classical", "Modified" };
const Teuchos::StringToIntMap OrthogonalizationMap(Orthogonalization_name,numOrthogonalizationValues,OrthogonalizationValues);
const EOrthogonalization Orthogonalization_default = ORTHOG_CLASSICAL;
// SizeOfKrylovSubspace
const int SizeOfKrylovSubspace_default = 300;
// ConvergenceTest
const int numConvergenceTestValues = 5;
enum EConvergenceTest
{ CONV_TEST_R0, CONV_TEST_RHS, CONV_TEST_ANORM, CONV_TEST_NOSCALING, CONV_TEST_SOL };
const char* ConvergenceTestValues[numConvergenceTestValues] =
{ "r0", "rhs", "Anorm", "no scaling", "sol" };
const Teuchos::StringToIntMap
ConvergenceTestMap(ConvergenceTest_name,numConvergenceTestValues,ConvergenceTestValues);
const EConvergenceTest ConvergenceTest_default = CONV_TEST_R0; // Same as RHS when X=0
// IllConditioningThreshold
const double IllConditioningThreshold_default = 1e+11;
// OutputFrequency
const int OutputFrequency_default = 0; // Aztec only prints to stdout which is not acceptable!

} // namespace

void setAztecOOParameters(
  Teuchos::ParameterList  *pl
  ,AztecOO                *solver
  )
{
  // Aztec Solver
  const std::string aztecSolverName = pl->get(AztecSolver_name,AztecSolverValues[AztecSolver_default]);
  switch(AztecSolverMap.get<EAztecSolver>(aztecSolverName,pl->name()+"->"+AztecSolver_name)) {
    case AZTEC_SOLVER_CG:
      solver->SetAztecOption(AZ_solver,AZ_cg);
      break;
    case AZTEC_SOLVER_GMRES:
      solver->SetAztecOption(AZ_solver,AZ_gmres);
      break;
    case AZTEC_SOLVER_CGS:
      solver->SetAztecOption(AZ_solver,AZ_cgs);
      break;
    case AZTEC_SOLVER_TFQMR:
      solver->SetAztecOption(AZ_solver,AZ_tfqmr);
      break;
    case AZTEC_SOLVER_BICGSTAB:
      solver->SetAztecOption(AZ_solver,AZ_bicgstab);
      break;
    case AZTEC_SOLVER_LU:
      solver->SetAztecOption(AZ_solver,AZ_lu);
      break;
    default:
      TEST_FOR_EXCEPT(true); // Should never get here!
  }
  // Aztec Preconditioner
  const::string aztecPrecName = pl->get(AztecPreconditioner_name,AztecPreconditionerValues[AztecPreconditioner_default]);
  switch(AztecPreconditionerMap.get<EAztecPreconditioner>(aztecPrecName,pl->name()+"->"+AztecPreconditioner_name)) {
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
  std::string rcmReorderingName = pl->get(RCMReordering_name,RCMReorderingValues[RCMReordering_default]);
  switch(RCMReorderingMap.get<ERCMReordering>(rcmReorderingName,pl->name()+"->"+RCMReordering_name)) {
    case RCM_REORDERING_ENABLED:
      solver->SetAztecOption(AZ_reorder,1);
      break;
    case RCM_REORDERING_DISABLED:
      solver->SetAztecOption(AZ_reorder,0);
      break;
    default:
      TEST_FOR_EXCEPT(true); // Should never get here!
  }
  // Gram-Schmidt orthogonalization procedure (GMRES only)
  std::string orthogName = pl->get(Orthogonalization_name,OrthogonalizationValues[Orthogonalization_default]);
  switch(OrthogonalizationMap.get<EOrthogonalization>(orthogName,pl->name()+"->"+Orthogonalization_name)) {
    case ORTHOG_CLASSICAL:
      solver->SetAztecOption(AZ_orthog,AZ_classic);
      break;
    case ORTHOG_MODIFIED:
      solver->SetAztecOption(AZ_orthog,AZ_modified);
      break;
    default:
      TEST_FOR_EXCEPT(true); // Should never get here!
  }
  // Size of the krylov subspace
  solver->SetAztecOption(AZ_kspace,pl->get(SizeOfKrylovSubspace_name,SizeOfKrylovSubspace_default));
  // Convergence criteria to use in the linear solver
  std::string convTestName = pl->get(ConvergenceTest_name,ConvergenceTestValues[ConvergenceTest_default]);
  switch(ConvergenceTestMap.get<EConvergenceTest>(convTestName,pl->name()+"->"+ConvergenceTest_name)) {
    case CONV_TEST_R0:
      solver->SetAztecOption(AZ_conv,AZ_r0);
      break;
    case CONV_TEST_RHS:
      solver->SetAztecOption(AZ_conv,AZ_rhs);
      break;
    case CONV_TEST_ANORM:
      solver->SetAztecOption(AZ_conv,AZ_Anorm);
      break;
    case CONV_TEST_NOSCALING:
      solver->SetAztecOption(AZ_conv,AZ_noscaled);
      break;
    case CONV_TEST_SOL:
      solver->SetAztecOption(AZ_conv,AZ_sol);
      break;
    default:
      TEST_FOR_EXCEPT(true); // Should never get here!
  }
  // Set the ill-conditioning threshold for the upper hessenberg matrix
  solver->SetAztecParam(AZ_ill_cond_thresh,pl->get(IllConditioningThreshold_name,IllConditioningThreshold_default));
  // Frequency of linear solve residual output
  solver->SetAztecOption(AZ_output,pl->get(OutputFrequency_name,OutputFrequency_default));
}

Teuchos::RefCountPtr<const Teuchos::ParameterList>
getValidAztecOOParameters()
{
  Teuchos::RefCountPtr<Teuchos::ParameterList>
    pl = Teuchos::rcp(new Teuchos::ParameterList());
  pl->set(AztecSolver_name,AztecSolverValues[AztecSolver_default]);
  pl->set(AztecPreconditioner_name,AztecPreconditionerValues[AztecPreconditioner_default]);
  pl->set(Overlap_name,Overlap_default);
  pl->set(GraphFill_name,GraphFill_default);
  pl->set(DropTolerance_name,DropTolerance_default);
  pl->set(FillFactor_name,FillFactor_default);
  pl->set(Steps_name,Steps_default);
  pl->set(PolynomialOrder_name,PolynomialOrder_default);
  pl->set(RCMReordering_name,RCMReorderingValues[RCMReordering_default]);
  pl->set(Orthogonalization_name,OrthogonalizationValues[Orthogonalization_default]);
  pl->set(SizeOfKrylovSubspace_name,SizeOfKrylovSubspace_default);
  pl->set(ConvergenceTest_name,ConvergenceTestValues[ConvergenceTest_default]);
  pl->set(IllConditioningThreshold_name,IllConditioningThreshold_default);
  pl->set(OutputFrequency_name,OutputFrequency_default);
  return pl;
}
