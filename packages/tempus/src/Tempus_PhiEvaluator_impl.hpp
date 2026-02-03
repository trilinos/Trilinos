//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_PhiEvaluator_impl_hpp
#define Tempus_PhiEvaluator_impl_hpp

#include "Tempus_PhiEvaluator.hpp"
#include "Tempus_PhiEvaluator_decl.hpp"
#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"

//#include "Thyra_DefaultInverseLinearOp_decl.hpp"
//#include "Thyra_DefaultMultipliedLinearOp_decl.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_SolveSupportTypes.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps_decl.hpp"
#include "Thyra_OperatorVectorTypes.hpp"
#include "Thyra_ProductVectorBase.hpp"

#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DefaultZeroLinearOp_decl.hpp"
#include "Thyra_DefaultZeroLinearOp_def.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"

#include "Thyra_DefaultScaledAdjointLinearOp.hpp"

// For printing
    #include "Teuchos_FancyOStream.hpp"
    #include "Teuchos_VerbosityLevel.hpp"


namespace Tempus {

template <class Scalar>
PhiEvaluator<Scalar>::PhiEvaluator()
  : name_("Phi Evaluator"), isInitialized_(false)
{
}

template <class Scalar>
PhiEvaluator<Scalar>::PhiEvaluator(std::string name)
  : isInitialized_(false)
{
  setName(name);
}

template <class Scalar>
void PhiEvaluator<Scalar>::copy(
    Teuchos::RCP<const PhiEvaluator<Scalar> > phi)
{
  this->setName(phi->getName());
}

template <class Scalar>
std::string PhiEvaluator<Scalar>::description() const
{
  return ("Tempus::PhiEvaluator - '" + name_ + "'");
}

template <class Scalar>
void PhiEvaluator<Scalar>::describe(
    Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const
{
  //TODO
  auto l_out = Teuchos::fancyOStream(out.getOStream());
  Teuchos::OSTab ostab(*l_out, 2, this->description());
  l_out->setOutputToRootOnly(0);

  *l_out << "\n--- " << this->description() << " ---" << std::endl;

  if ((Teuchos::as<int>(verbLevel) ==
       Teuchos::as<int>(Teuchos::VERB_DEFAULT)) ||
      (Teuchos::as<int>(verbLevel) >= Teuchos::as<int>(Teuchos::VERB_LOW))) {
    //*l_out << "  abc     = " << ... <<
    // std::endl;
  }

  //if (Teuchos::as<int>(verbLevel) >= Teuchos::as<int>(Teuchos::VERB_MEDIUM)) {
  //
  //}
  *l_out << std::string(this->description().length() + 8, '-') << std::endl;
}

template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
PhiEvaluator<Scalar>::getValidParameters() const
{
  return this->getValidParametersBasic();
}

template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
PhiEvaluator<Scalar>::getValidParametersBasic() const
{
  //TODO
  Teuchos::RCP<Teuchos::ParameterList> pl =
      Teuchos::parameterList("Phi Evaluator");

  pl->setName(getName());

  pl->set(
      "PhiEvaluator Type", "PFD",
      "Method to approximate the phi-function evaluation.");

  pl->set(
      "Lump Mass", false,
      "'Lump Mass' switch on lumping of the mass matrix.  "
      "'true' - will lump the mass matrix for PhiEvaluators that support this feature.  "
      "'false' - will switch off mass lumping.");

  //pl->set("?", *member_->getNonconstParameterList());

  return pl;
}

template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
PhiEvaluator<Scalar>::getNonconstParameterList()
{
  return Teuchos::rcp_const_cast<Teuchos::ParameterList>(getValidParameters());
}

template <class Scalar>
void PhiEvaluator<Scalar>::initialize() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      appModel_ == Teuchos::null, std::logic_error,
      "Error - PhiEvaluator::initialize() Model not set!\n");

  isInitialized_ = true;  // Only place where this is set to true!
}

template <class Scalar>
void PhiEvaluator<Scalar>::checkInitialized()
{
  if (!this->isInitialized()) {
    this->describe(*(this->getOStream()), Teuchos::VERB_MEDIUM);
    TEUCHOS_TEST_FOR_EXCEPTION(
        !this->isInitialized(), std::logic_error,
        "Error - " << this->description() << " is not initialized!");
  }
}

template <class Scalar>
void PhiEvaluator<Scalar>::setModel(const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > appModel)
{
  appModel_ = appModel;

  // TODO: read this in from the xml file
  const bool lumpmass = true;
  phiLinSolv_ = Teuchos::rcp(new PhiLinearSolver<Scalar>(appModel_, lumpmass));
}


/*
 * PhiLinearSolver methods
 */

template <class Scalar>
void PhiLinearSolver<Scalar>::computeMassMatrix(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs)
// void PhiLinearSolver<Scalar>::computeMassMatrix(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs)
{
  
  typedef Thyra::ModelEvaluatorBase MEB;

  // MEB::InArgs<Scalar> inArgs = appModel_->getNominalValues();

  // first allocate space for the mass matrix
  fullMassMatrix_ = appModel_->create_W_op();

  // request only the mass matrix from the physics
  // Model evaluator builds: alpha*M*u_dot + beta*F(u) = 0
  MEB::InArgs<Scalar> inArgs_new  = appModel_->createInArgs();
  inArgs_new.setArgs(inArgs);
  inArgs_new.set_x_dot(inArgs.get_x_dot()); //TODO: why?
  inArgs_new.set_alpha(1.0);
  inArgs_new.set_beta(0.0);

  // TODO:
  // set the one time beta to ensure dirichlet conditions
  // are correctly included in the mass matrix: do it for
  // both epetra and Tpetra.
  //if(panzerModel_!=Teuchos::null)
  //  panzerModel_->setOneTimeDirichletBeta(1.0);
  //else {
    // assuming the underlying model is a delegator, walk through
    // the decerator hierarchy until you find a panzer::ME or panzer::EpetraME.
    // If you don't find one, then throw because you are in a load of trouble anyway!
  //  setOneTimeDirichletBeta(1.0,*this->getUnderlyingModel());
  //}

  // set only the mass matrix
  MEB::OutArgs<Scalar> outArgs = appModel_->createOutArgs();
  outArgs.set_W_op(fullMassMatrix_);

  // this will fill the mass matrix operator
  appModel_->evalModel(inArgs_new, outArgs);

  //TODO:
  //Teuchos::RCP<const Epetra_CrsMatrix> crsMat = Teuchos::rcp_dynamic_cast<const Epetra_CrsMatrix>(Thyra::get_Epetra_Operator(*fullMassMatrix));
  //EpetraExt::RowMatrixToMatrixMarketFile("fullMassMatrix_mat.mm",*crsMat);

  if(!lumpMass_) {
    invMassMatrix_ = Thyra::inverse<Scalar>(*appModel_->get_W_factory(),fullMassMatrix_);
  }
  else {
    // build lumped mass matrix (assumes all positive mass entries, does a simple sum)
    Teuchos::RCP<Thyra::VectorBase<Scalar> > ones = Thyra::createMember(*fullMassMatrix_->domain());
    Thyra::assign(ones.ptr(),1.0);

    lumpedMassDiagonal_ = Thyra::createMember(*fullMassMatrix_->range());
    Teuchos::RCP<Thyra::VectorBase<Scalar> > invLumpMass = Thyra::createMember(*fullMassMatrix_->range());
    Thyra::apply(*fullMassMatrix_, Thyra::NOTRANS, *ones, lumpedMassDiagonal_.ptr());

    //TODO:
    //Teuchos::RCP<const Epetra_Vector> mv = Teuchos::rcp_dynamic_cast<const Epetra_Vector>(Thyra::get_Epetra_Vector(crsMat->RangeMap(), invLumpMass));
    //EpetraExt::VectorToMatrixMarketFile("fullMassMatrix_v.mm", *mv);

    Thyra::reciprocal(*lumpedMassDiagonal_, invLumpMass.ptr());

    lumpMassMatrix_ = Thyra::diagonal(lumpedMassDiagonal_);
    invMassMatrix_ = Thyra::diagonal(invLumpMass);
  }
}

template <class Scalar>
Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> PhiLinearSolver<Scalar>::buildATilde(
    const double dt)   // time step
{
  // Combine linear operators M_inv and J and multiply by -dt (minus is for implicit to explicit conversion)
  Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> A = Thyra::scale(Teuchos::as<Scalar>(-dt), Thyra::multiply<Scalar>(invMassMatrix_, jacobianMatrix_));

  // Spaces
  auto V = A->domain();   // dim N, and also A->range()
  auto W = KMatrix_->domain();   // dim p, and also K->range()
  
  // Zero operator: V -> W  (for the (2,1) block)
  auto Z_VW = Thyra::zero<Scalar>(W, V);

  // Build block operator:
  // [ A  b ]
  // [ 0  K ]
  Atilde_ = Thyra::block2x2<Scalar>(
      A,            // (1,1): V->V
      bMatrix_,     // (1,2): W->V
      Z_VW,         // (2,1): V->W
      KMatrix_      // (2,2): W->W
  );

  return Atilde_;
}

template <class Scalar>
void PhiLinearSolver<Scalar>::buildK(const Thyra::Ordinal p)
{
  // Space of dimension p
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>> V =
      Thyra::defaultSpmdVectorSpace<Scalar>(p);

  // Create a p-column multivector: rows = p, cols = p
  Teuchos::RCP<Thyra::MultiVectorBase<Scalar>> K_mv = Thyra::createMembers(V, p);

  // Initialize to zero
  Thyra::assign(K_mv.ptr(), Scalar(0));

  // Fill superdiagonal: K(i, i+1) = 1
  for (Thyra::Ordinal j = 1; j < p; ++j)
  {
      // Column j, row j-1
      auto col_j = K_mv->col(j); 
      Thyra::set_ele(j - 1, Scalar(1), col_j.ptr());
  }

  // Wrap as LinearOp
  KMatrix_ = K_mv;
}

template <class Scalar>
void PhiLinearSolver<Scalar>::buildb(const Thyra::Ordinal p,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& xDot)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      p > 2,
      std::invalid_argument,
      "buildb: p must be 2. Higher order EPI is not yet supported.");

  TEUCHOS_TEST_FOR_EXCEPTION(
      p < 2,
      std::invalid_argument,
      "buildb: EPI order must be 2 or higher.");

  // N-dimensional space: use A's range (rows of an NxN operator)
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>> V_N = invMassMatrix_->range();

  const Thyra::Ordinal N = V_N->dim();

  // Create an N x p multivector (N rows, p columns)
  Teuchos::RCP<Thyra::MultiVectorBase<Scalar>> b_Np = Thyra::createMembers(V_N, p);

  // Initialize to zero
  Thyra::assign(b_Np.ptr(), Scalar(0));

  // Fill the 2nd column with xDot (Frhs) TODO: This needs to be updated for higher order support
  auto col1 = b_Np->col(1);
  Thyra::assign(col1.ptr(), *xDot);

  // Store b
  bMatrix_ = b_Np;
}

template <class Scalar>
Teuchos::RCP<Thyra::VectorBase<Scalar>> PhiLinearSolver<Scalar>::buildv(const Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>> space)
{
  // v must be in the domain of Atilde_
  const Thyra::Ordinal dim = space->dim();

  // Create v and initialize to zero
  v_ = Thyra::createMember(space);
  Thyra::assign(v_.ptr(), Scalar(0));

  // Get the last index
  const Thyra::Ordinal g_last = dim - 1;

  // If this is an distributed space, set only on the owning rank
  if (auto spmdSpace = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<Scalar>>(space))
  {
    const Thyra::Ordinal localOffset = spmdSpace->localOffset();
    const Thyra::Ordinal localSubDim = spmdSpace->localSubDim();

    if (g_last >= localOffset && g_last < localOffset + localSubDim)
      Thyra::set_ele(g_last, Scalar(1), v_.ptr());

  }
  else
  {
    // Fallback for non-SPMD spaces
    Thyra::set_ele(g_last, Scalar(1), v_.ptr());
  }

  return v_;
}

// template <class Scalar>
// Teuchos::RCP<const Thyra::VectorBase<Scalar>> PhiLinearSolver<Scalar>::matrixExponential(const int expansionOrder)
// {
//   TEUCHOS_TEST_FOR_EXCEPTION(
//       expansionOrder < 0,
//       std::invalid_argument,
//       "matrixExponential: expansionOrder must be nonnegative");

//   // exp(A) * v is in range(A)
//   const auto rangeSpace = Atilde_->range();

//   // Create tmp vector to hold result
//   auto matExpTemp = Thyra::createMember(rangeSpace);

//   Thyra::assign(matExpTemp.ptr(), Scalar(0));

//   // Identity * v = v
//   Teuchos::RCP<Thyra::VectorBase<Scalar>> term = Thyra::createMember(rangeSpace);
//   Thyra::assign(term.ptr(), *v_);

//   // matExpTemp += term / 0!
//   Thyra::Vp_V(matExpTemp.ptr(), *term);

//   // Iteratively compute term = A * term (A^k v) and accumulate term/k!
//   Scalar invFact = Scalar(1); // 1/k! updated each step
//   for (int k = 1; k <= expansionOrder; ++k)
//   {
//     // term <- A * term
//     Teuchos::RCP<Thyra::VectorBase<Scalar>> next = Thyra::createMember(rangeSpace);
//     Thyra::apply(*Atilde_, Thyra::NOTRANS, *term, next.ptr());
//     term = next;

//     invFact /= Scalar(k);

//     // multiply with inverse factorial
//     Thyra::Vp_StV(matExpTemp.ptr(), invFact, *term);
//   }

//   matExp_v_ = matExpTemp; // This is required to wrap multivector as linearop

//   return matExp_v_;
// }


template <class Scalar>
void PhiLinearSolver<Scalar>::applyMass(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> Mf,
                                        const Teuchos::RCP<const Thyra::VectorBase<Scalar>> f) const
{
  // apply the mass matrix
  if (f != Teuchos::null && Mf != Teuchos::null) {
    if (!lumpMass_)
      Thyra::apply(*fullMassMatrix_, Thyra::NOTRANS, *f, Mf.ptr());
    else
      Thyra::apply(*lumpMassMatrix_, Thyra::NOTRANS, *f, Mf.ptr());
  }
}

template <class Scalar>
void PhiLinearSolver<Scalar>::solveMass(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> f,
                                        const Teuchos::RCP<const Thyra::VectorBase<Scalar>> Mf) const
{
  // invert the mass matrix
  if (Mf != Teuchos::null && f != Teuchos::null) {
    Thyra::apply(*invMassMatrix_, Thyra::NOTRANS, *Mf, f.ptr());
  }
}

template <class Scalar>
void PhiLinearSolver<Scalar>::computeJacobian(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs)
{
  typedef Thyra::ModelEvaluatorBase MEB;

  // MEB::InArgs<Scalar> inArgs = appModel_->getNominalValues();

  // first allocate space for the Jacobian matrix
  jacobianMatrix_ = appModel_->create_W_op();

  // request only the Jacobian matrix from the physics
  // Model evaluator builds: alpha*u_dot + beta*F(u) = 0
  MEB::InArgs<Scalar> inArgs_new  = appModel_->createInArgs();
  inArgs_new.setArgs(inArgs);
  inArgs_new.set_x_dot(inArgs.get_x_dot());
  inArgs_new.set_alpha(0.0);
  inArgs_new.set_beta(1.0);

  // set only the Jacobian matrix
  MEB::OutArgs<Scalar> outArgs = appModel_->createOutArgs();
  outArgs.set_W_op(jacobianMatrix_);

  // this will fill the Jacobian matrix operator
  appModel_->evalModel(inArgs_new, outArgs);
}

template <class Scalar>
void PhiLinearSolver<Scalar>::applyJacobian(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> Jf,
                                            const Teuchos::RCP<const Thyra::VectorBase<Scalar>> f) const
{
  // apply the Jacobian matrix
  if (f != Teuchos::null && Jf != Teuchos::null) {
    Thyra::apply(*jacobianMatrix_, Thyra::NOTRANS, *f, Jf.ptr());
  }
}

template <class Scalar>
Thyra::SolveStatus<Scalar> PhiLinearSolver<Scalar>::solveMpJ(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
							     const Teuchos::Ptr<Thyra::VectorBase<Scalar>> iMf,
							     const Teuchos::RCP<const Thyra::VectorBase<Scalar>> Mf,
							     Scalar alpha, Scalar beta) const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  MEB::OutArgs<Scalar> outArgs = appModel_->createOutArgs();

  // first allocate space for the jacobian
  Teuchos::RCP<Thyra::LinearOpBase<Scalar>> MpJ = appModel_->create_W_op();
  Teuchos::RCP<Thyra::PreconditionerBase<Scalar>> MpJ_p = Teuchos::null;
  if (outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_W_prec)) {
    MpJ_p = appModel_->create_W_prec();
  }
  //TODO: support all types of preconditoner (compare NOX_Thyra_Group.C)

  MEB::InArgs<Scalar> inArgs_new = appModel_->createInArgs();
  inArgs_new.setArgs(inArgs);
  inArgs_new.set_x_dot(inArgs.get_x_dot());
  inArgs_new.set_alpha(alpha);
  inArgs_new.set_beta(beta);

  // set only the mass plus jacobian (MpJ) matrix
  outArgs.set_W_op(MpJ);
  if (MpJ_p != Teuchos::null){
    outArgs.set_W_prec(MpJ_p);
  }

  // this will fill the MpJ operator
  appModel_->evalModel(inArgs_new, outArgs);

  // TODO: const-cast why?
  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar>> const_lowsFactory = appModel_->get_W_factory();
  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar>> lowsFactory =
    Teuchos::rcp_const_cast<Thyra::LinearOpWithSolveFactoryBase<Scalar>>(const_lowsFactory);
      
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar>> LOWSB = Teuchos::null;
  if (MpJ_p == Teuchos::null){
    // without preconditioner
    Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> const_MpJ = Teuchos::rcpFromRef(*MpJ);
    LOWSB = Thyra::linearOpWithSolve(*lowsFactory, const_MpJ);
  }
  else {
    // with preconditioner
    LOWSB = lowsFactory->createOp();
    Thyra::initializePreconditionedOp<Scalar>(*lowsFactory, MpJ, MpJ_p, LOWSB.ptr());
  }

  assign(iMf.ptr(), ST::zero()); //TODO: needed?

  // Create solve criteria
  Thyra::SolveCriteria<Scalar> solveCriteria;
  solveCriteria.requestedTol = 1e-6;//p.get("Tolerance", 1.0e-6);

  //std::string numer_measure = p.get("Solve Measure Numerator", "Norm Residual");
  //std::string denom_measure = p.get("Solve Measure Denominator", "Norm Initial Residual");

  //if (name == "None")
  //  return ::Thyra::SOLVE_MEASURE_ONE;
  //else if (name == "Norm Residual")
  //  return ::Thyra::SOLVE_MEASURE_NORM_RESIDUAL;
  //else if (name == "Norm Solution")
  //  return ::Thyra::SOLVE_MEASURE_NORM_SOLUTION;
  //else if (name == "Norm Initial Residual")
  //  return ::Thyra::SOLVE_MEASURE_NORM_INIT_RESIDUAL;
  //else if (name == "Norm RHS")
  //  return ::Thyra::SOLVE_MEASURE_NORM_RHS;

  solveCriteria.solveMeasureType =
    Thyra::SolveMeasureType(Thyra::SOLVE_MEASURE_NORM_RESIDUAL, Thyra::SOLVE_MEASURE_NORM_INIT_RESIDUAL);

  // compute the solution to (MpJ)\Mf and write it to iMf
  Thyra::SolveStatus<Scalar> sStatus = LOWSB->solve(Thyra::NOTRANS, *Mf, iMf.ptr(), Teuchos::constPtr(solveCriteria));

  return sStatus;
}


}  // namespace Tempus
#endif  // Tempus_PhiEvaluator_impl_hpp

/**
 * TODO: code for working with preconditioners
      using Teuchos::tab; using Teuchos::rcpFromPtr;
      typedef Teuchos::ScalarTraits<Scalar> ST;
      Teuchos::OSTab tab2(out);
      out << "\nShowing resuse of the preconditioner ...\n";
      Teuchos::RCP<Thyra::LinearOpBase<Scalar> > A = rcpFromPtr(A_inout);
      // Create the initial preconditioner for the input forward operator
      Teuchos::RCP<Thyra::PreconditionerBase<Scalar> > P =
      precFactory.createPrec();
      Thyra::initializePrec<Scalar>(precFactory, A, P.ptr());
      // Create the invertible LOWS object given the preconditioner
      Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> > invertibleA =
      lowsFactory.createOp();
      Thyra::initializePreconditionedOp<Scalar>(lowsFactory, A, P, invertibleA.ptr());
      // Solve the first linear system
      assign(x1, ST::zero());
      Thyra::SolveStatus<Scalar> status1 = Thyra::solve<Scalar>(*invertibleA,
      Thyra::NOTRANS, b1, x1);
      out << "\nSolve status:\n" << status1;
      // Change the forward linear operator without changing the preconditioner
      opChanger.changeOp(A.ptr());
      // Warning! After the above change the integrity of the preconditioner
      // linear operators in P is undefined. For some implementations of the
      // preconditioner, its behavior will remain unchanged (e.g. ILU) which in
      // other cases the behavior will change but the preconditioner will still
      // work (e.g. Jacobi). However, there may be valid implementations where
      // the preconditioner will simply break if the forward operator that it is
      // based on breaks.
      //
      // Reinitialize the LOWS object given the updated forward operator A and the
      // old preconditioner P.
      Thyra::initializePreconditionedOp<Scalar>(lowsFactory, A, P, invertibleA.ptr());
      // Solve the second linear system
      assign(x2, ST::zero());
      Thyra::SolveStatus<Scalar>status2 = Thyra::solve<Scalar>(*invertibleA,
      Thyra::NOTRANS, b2, x2);
      out << "\nSolve status:\n" << status2;
      } // end externalPreconditionerReuseWithSolves
 */



/**
 * TODO: code from Panzer explicit modelevaluator
 *
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"

#include "Panzer_ExplicitModelEvaluator.hpp"

#include "PanzerDiscFE_config.hpp"

namespace panzer {

template<typename Scalar>
ExplicitModelEvaluator<Scalar>::
ExplicitModelEvaluator(const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > & model,
                       bool constantMassMatrix,
                       bool useLumpedMass,
                       bool applyMassInverse)
   : Thyra::ModelEvaluatorDelegatorBase<Scalar>(model)
   , constantMassMatrix_(constantMassMatrix)
   , massLumping_(useLumpedMass)
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  this->applyMassInverse_ = applyMassInverse;

  // extract a panzer::ModelEvaluator if appropriate
  panzerModel_ = rcp_dynamic_cast<panzer::ModelEvaluator<Scalar> >(model);

  // note at this point its possible that panzerModel_ = panzerEpetraModel_ = Teuchos::null

  buildArgsPrototypes();
}

template<typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar> ExplicitModelEvaluator<Scalar>::
getNominalValues() const
{
  typedef Thyra::ModelEvaluatorBase MEB;

  MEB::InArgs<Scalar> nomVals = createInArgs();
  nomVals.setArgs(this->getUnderlyingModel()->getNominalValues(),true);

  return nomVals;
}

template<typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar> ExplicitModelEvaluator<Scalar>::
createInArgs() const
{
  return prototypeInArgs_;
}

template<typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar> ExplicitModelEvaluator<Scalar>::
createOutArgs() const
{
  return prototypeOutArgs_;
}

template<typename Scalar>
Teuchos::RCP<panzer::ModelEvaluator<Scalar> > ExplicitModelEvaluator<Scalar>::
getPanzerUnderlyingModel()
{
  return Teuchos::rcp_dynamic_cast<panzer::ModelEvaluator<Scalar> >(this->getNonconstUnderlyingModel());
}

template<typename Scalar>
void ExplicitModelEvaluator<Scalar>::
evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
              const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  using Teuchos::RCP;
  RCP<const Thyra::ModelEvaluator<Scalar> > under_me = this->getUnderlyingModel();

  MEB::InArgs<Scalar> under_inArgs = under_me->createInArgs();
  under_inArgs.setArgs(inArgs);

  // read in the supplied time derivative in case it is needed by the explicit model (e.g. for stabilization) and make sure alpha is set to zero
  under_inArgs.set_x_dot(inArgs.get_x_dot());
  under_inArgs.set_alpha(0.0);

  Teuchos::RCP<Thyra::VectorBase<Scalar> > f = outArgs.get_f();

  MEB::OutArgs<Scalar> under_outArgs = under_me->createOutArgs();
  under_outArgs.setArgs(outArgs);
  if(f!=Teuchos::null) {
    // build a scrap vector that will contain the weak residual
    if(scrap_f_==Teuchos::null)
      scrap_f_ = Thyra::createMember(*under_me->get_f_space());

    Thyra::assign(scrap_f_.ptr(),0.0);
    under_outArgs.set_f(scrap_f_);
  }

  under_me->evalModel(under_inArgs,under_outArgs);

  // build the mass matrix
  if(invMassMatrix_==Teuchos::null || constantMassMatrix_==false)
    buildInverseMassMatrix(inArgs);

  if(f!=Teuchos::null)
      Thyra::Vt_S(scrap_f_.ptr(),-1.0);

  // invert the mass matrix
  if(f!=Teuchos::null && this->applyMassInverse_) {
    Thyra::apply(*invMassMatrix_,Thyra::NOTRANS,*scrap_f_,f.ptr());
  }
  else if(f!=Teuchos::null){
    Thyra::V_V(f.ptr(),*scrap_f_);
  }
}

template<typename Scalar>
void ExplicitModelEvaluator<Scalar>::
buildInverseMassMatrix(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs) const
{
}

template<typename Scalar>
void ExplicitModelEvaluator<Scalar>::
buildArgsPrototypes()
{
  typedef Thyra::ModelEvaluatorBase MEB;

  MEB::InArgsSetup<Scalar> inArgs(this->getUnderlyingModel()->createInArgs());
  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(MEB::IN_ARG_alpha,true);
  inArgs.setSupports(MEB::IN_ARG_beta,true);
  inArgs.setSupports(MEB::IN_ARG_x_dot,true);
  prototypeInArgs_ = inArgs;

  MEB::OutArgsSetup<Scalar> outArgs(this->getUnderlyingModel()->createOutArgs());
  outArgs.setModelEvalDescription(this->description());
  outArgs.setSupports(MEB::OUT_ARG_W,false);
  outArgs.setSupports(MEB::OUT_ARG_W_op,false);
  prototypeOutArgs_ = outArgs;
}

template<typename Scalar>
void ExplicitModelEvaluator<Scalar>::
setOneTimeDirichletBeta(double beta,const Thyra::ModelEvaluator<Scalar> & me) const
{
  using Teuchos::Ptr;
  using Teuchos::ptrFromRef;
  using Teuchos::ptr_dynamic_cast;

  // try to extract base classes that support setOneTimeDirichletBeta
  Ptr<const panzer::ModelEvaluator<Scalar> > panzerModel = ptr_dynamic_cast<const panzer::ModelEvaluator<Scalar> >(ptrFromRef(me));
  if(panzerModel!=Teuchos::null) {
    panzerModel->setOneTimeDirichletBeta(beta);
    return;
  }

  // if you get here then the ME is not a panzer::ME or panzer::EpetraME, check
  // to see if its a delegator

  Ptr<const Thyra::ModelEvaluatorDelegatorBase<Scalar> > delegator
      = ptr_dynamic_cast<const Thyra::ModelEvaluatorDelegatorBase<Scalar> >(ptrFromRef(me));
  if(delegator!=Teuchos::null) {
    setOneTimeDirichletBeta(beta,*delegator->getUnderlyingModel());
    return;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                             "panzer::ExplicitModelEvaluator::setOneTimeDirichletBeta can't find a panzer::ME or a panzer::EpetraME. "
                             "The deepest model is also not a delegator. Thus the recursion failed and an exception was generated.");
}

} // end namespace panzer

#endif
**/
