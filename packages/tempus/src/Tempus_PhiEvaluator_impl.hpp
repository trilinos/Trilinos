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
#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Tempus_config.hpp"

#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"

#include "Thyra_DefaultInverseLinearOp.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_SolveSupportTypes.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_OperatorVectorTypes.hpp"
#include "Thyra_ProductVectorBase.hpp"

#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"

#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultAddedLinearOp_def.hpp"

// For printing
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_VerbosityLevel.hpp"

// Anasazi eigensolver (Block Krylov-Schur) for spectrum bound estimation
#ifdef TEMPUS_ENABLE_ANASAZI
#include "AnasaziThyraAdapter.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#endif


namespace Tempus {

template <class Scalar>
PhiEvaluator<Scalar>::PhiEvaluator()
  : PhiEvaluator<Scalar>("Phi Evaluator")
{
}

template <class Scalar>
PhiEvaluator<Scalar>::PhiEvaluator(std::string name)
  : isInitialized_(false),
    lumpMassMatrix_(false),
    constantMassMatrix_(false),
    useAtildeForSingleRHS_(true) // TODO: make this configurable
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
  auto l_out = Teuchos::fancyOStream(out.getOStream());
  Teuchos::OSTab ostab(*l_out, 2, this->description());
  l_out->setOutputToRootOnly(0);

  *l_out << "\n--- " << this->description() << " ---" << std::endl;

  if ((Teuchos::as<int>(verbLevel) ==
       Teuchos::as<int>(Teuchos::VERB_DEFAULT)) ||
      (Teuchos::as<int>(verbLevel) >= Teuchos::as<int>(Teuchos::VERB_LOW))) {
    *l_out << "  lumpMassMatrix_        = " << lumpMassMatrix_ << std::endl;
    *l_out << "  constantMassMatrix_    = " << constantMassMatrix_ << std::endl;
    *l_out << "  phiLinSolv_            = " << phiLinSolv_ << std::endl;
    *l_out << "  useAtildeForSingleRHS_ = " << useAtildeForSingleRHS_ << std::endl;
  }

  //if (Teuchos::as<int>(verbLevel) >= Teuchos::as<int>(Teuchos::VERB_MEDIUM)) {
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
  Teuchos::RCP<Teuchos::ParameterList> pl =
      Teuchos::parameterList("Phi Evaluator");

  pl->setName(getName());

  pl->set(
      "PhiEvaluator Type", "PFD",
      "Method to approximate the phi-function evaluation.");

  pl->set(
      "Lump Mass Matrix", false,
      "'Lump Mass Matrix' switch on lumping of the mass matrix."
      "'true' - will lump the mass matrix for PhiEvaluators that support this feature."
      "'false' - will switch off mass lumping.");

  pl->set(
      "Constant Mass Matrix", true,
      "'Constant Mass Matrix' caches the mass matrix, which is assumed constant in time."
      "'true' - will cache the mass matrix."
      "'false' - will recompute the mass matrix every time-step.");

  pl->sublist("Eigensolver");
  pl->sublist("Eigensolver").set("Which", "LM");
  pl->sublist("Eigensolver").set("Block Size",                     2);
  pl->sublist("Eigensolver").set("Num Blocks",                     16);
  pl->sublist("Eigensolver").set("Maximum Restarts",               50);
  pl->sublist("Eigensolver").set("Convergence Tolerance",          0.08);
  pl->sublist("Eigensolver").set("Relative Convergence Tolerance", true);
  pl->sublist("Eigensolver").set("Num Eigenvalues",                8,
      "Number of eigenvalues used to estimate the Jacobian spectrum bounds.");
  pl->sublist("Eigensolver").set("Dense Fallback Dimension",       64,
      "Systems with dimension <= this value use a dense LAPACK eigensolver "
      "instead of the iterative Block Krylov-Schur solver.");

  

  return pl;
}

template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
PhiEvaluator<Scalar>::getNonconstParameterList()
{
  return Teuchos::rcp_const_cast<Teuchos::ParameterList>(getValidParameters());
}

template <class Scalar>
void PhiEvaluator<Scalar>::initialize()
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      appModel_ == Teuchos::null, std::logic_error,
      "Error - PhiEvaluator::initialize() Model not set!\n");
  TEUCHOS_TEST_FOR_EXCEPTION(
      phiLinSolv_ == Teuchos::null, std::logic_error,
      "Error - PhiEvaluator::initialize() phiLinearSolver not set!\n");

  isInitialized_ = true;  // Only place where this is set to true!
}

template <class Scalar>
void PhiEvaluator<Scalar>::setPhiEvaluatorValues(
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  pl->validateParametersAndSetDefaults(*getValidParameters());

  setLumpMassMatrix(pl->get<bool>("Lump Mass Matrix", false));
  setConstantMassMatrix(pl->get<bool>("Constant Mass Matrix", false));
  setName(pl->name());

  // Cache a copy of the "Eigensolver" sublist so setModel() can forward it
  // to the PhiLinearSolver when it is constructed.
  eigensolverPL_ =
      Teuchos::rcp(new Teuchos::ParameterList(pl->sublist("Eigensolver")));
  if (phiLinSolv_ != Teuchos::null)
    phiLinSolv_->setEigensolverParams(eigensolverPL_);
}

template <class Scalar>
void PhiEvaluator<Scalar>::checkInitialized() const
{
  if (!this->isInitialized()) {
    this->describe(*(this->getOStream()), Teuchos::VERB_MEDIUM);
    TEUCHOS_TEST_FOR_EXCEPTION(
        !this->isInitialized(), std::logic_error,
        "Error - " << this->description() << " is not initialized!");
  }
}

template <class Scalar>
void PhiEvaluator<Scalar>::setModel(const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>> appModel)
{
  appModel_ = appModel;

  phiLinSolv_ = Teuchos::rcp(new PhiLinearSolver<Scalar>(appModel_, lumpMassMatrix_));

  // Forward eigensolver parameters if they were set before setModel was called
  if (eigensolverPL_ != Teuchos::null)
    phiLinSolv_->setEigensolverParams(eigensolverPL_);
}

template<class Scalar>
void PhiEvaluator<Scalar>::setLumpMassMatrix(bool lumpMassMatrix)
{
  lumpMassMatrix_ = lumpMassMatrix;
  if (this->phiLinSolv_ != Teuchos::null) {
    this->phiLinSolv_->setLumpMassMatrix(lumpMassMatrix_);
  }
}

template<class Scalar>
void PhiEvaluator<Scalar>::setConstantMassMatrix(bool constantMassMatrix)
{
  if (!constantMassMatrix && constantMassMatrix && this->phiLinSolv_ != Teuchos::null) {
    this->phiLinSolv_->clearMemory();
  }
  constantMassMatrix_ = constantMassMatrix;
}

template <class Scalar>
void PhiEvaluator<Scalar>::setLinearizationPoint(const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
                                                 const PhiInitialization& mode)
{
  // TODO: remove this copy, when the TakeStep goes out of scope, this pointer becomes invalid
  //       it is only used if we need inArgs after this method finished, i.e., need to assemble more matrices later.
  inArgs_lin_ = Teuchos::rcpFromRef(inArgs);

  // ensure that model and phiLinSolv are available
  checkInitialized();

  //if (mode == PhiInitialization::JACOBIAN_AND_MASS) {
    TEMPUS_FUNC_TIME_MONITOR("Tempus::PhiEvaluator::setLinearizationPoint (Mass and Jac assembly)");
    {
      // compute all required matrices at current linearization point
      this->phiLinSolv_->setLumpMassMatrix(this->lumpMassMatrix_);
      this->phiLinSolv_->computeMassMatrix(*inArgs_lin_);
      this->phiLinSolv_->computeJacobian(*inArgs_lin_);
    }
  //}
  //else {
  //  TEMPUS_FUNC_TIME_MONITOR("Tempus::PhiEvaluator::setLinearizationPoint (Mass assembly)");
  //  {
  //    // compute all required matrices at current linearization point
  //    this->phiLinSolv_->setLumpMassMatrix(this->lumpMassMatrix_);
  //    this->phiLinSolv_->computeMassMatrix(*inArgs_lin_);
  //  }
  //}
}

template<class Scalar>
Thyra::SolveStatus<Scalar>
PhiEvaluator<Scalar>::computePhi(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> x,
                                 const int phi_order, const Scalar cdt,
                                 const Teuchos::RCP<const Thyra::VectorBase<Scalar>> &Mrhs_b)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    phi_order < 0, std::logic_error,
    "Error - PhiEvaluator::computePhi() phi_order must be non-negative!\n");

  if (useAtildeForSingleRHS_ && phi_order > 0)
  {
    // Use linear combination extension formula and matrix exponential
    Teuchos::ArrayRCP<Teuchos::RCP<const Thyra::VectorBase<Scalar>>> Mrhs_B(phi_order+1);
    Mrhs_B[phi_order] = Mrhs_b;
    return this->computePhis(x, cdt, Mrhs_B());
  }
  else
  {
    // Call LinOpPhi directly

    // check that everything has been initialized properly
    checkInitialized();

    // Invert the mass matrix out of the right hand side and store in x
    Thyra::assign(x, 0.0);
    this->phiLinSolv_->solveMass(x, Mrhs_b);

    // Build linop with Jacobian and mass matrix
    const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> Lop = this->phiLinSolv_->buildL(cdt);

    // Compute phi-function in place
    Thyra::SolveStatus<Scalar> sStatus = this->computeLinOpPhi(phi_order, Lop, x, cdt);
    return sStatus;
  }
}

template <class Scalar>
Thyra::SolveStatus<Scalar>
PhiEvaluator<Scalar>::computePhis(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> x,
                                  const Scalar cdt,
                                  const Teuchos::ArrayView<const Teuchos::RCP<const Thyra::VectorBase<Scalar>>> &Mrhs_B)
{
  const Thyra::Ordinal max_phi_order = Mrhs_B.size() - 1;

  TEUCHOS_TEST_FOR_EXCEPTION(
      max_phi_order < 0,
      std::invalid_argument,
      "Error - PhiEvaluator::computePhis() list of rhs must have at least one entry.");

  // support max_phi_order == 0 by calling the non-extended solver
  if (max_phi_order == 0)
  {
    return computePhi(x, 0, cdt, Mrhs_B[0]);
  }

  // check that everything has been initialized properly
  checkInitialized();

  Teuchos::ArrayRCP<Teuchos::RCP<const Thyra::VectorBase<Scalar>>> rhs_B(max_phi_order+1);

  // Invert the mass matrix out of the right hand sides
  // TODO: This might be more efficient to do on the combined MultiVector that will be assembled in buildb
  //       However, if Mrhs_B is sparse, it may not.
  for (Thyra::Ordinal ii = 0; ii < max_phi_order+1; ii++)
  {
    if (Mrhs_B[ii] != Teuchos::null)
    {
      auto Mrhs_b = Mrhs_B[ii];
      Teuchos::RCP<Thyra::VectorBase<Scalar>> rhs_b = Mrhs_b->clone_v();
      Thyra::assign(rhs_b.ptr(), 0.0);
      this->phiLinSolv_->solveMass(rhs_b.ptr(), Mrhs_b);
      rhs_B[ii] = rhs_b;
    }
  }

  // Build extended matrix
  const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> ATilde
    = this->phiLinSolv_->buildATilde(cdt, rhs_B());

  // Build initial vector and compute matrix exponential in place
  auto v = this->phiLinSolv_->buildv(ATilde->domain(), rhs_B[0]);

  //Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  //Atilde->describe(*out, Teuchos::VERB_EXTREME);
  //v->describe(*out, Teuchos::VERB_EXTREME);

  Thyra::SolveStatus<Scalar> sStatus = this->computeLinOpPhi(0, ATilde, v.ptr(), cdt);

  // Get the first block of the multi-vector calculated from 2x2 multi-matrix
  auto v0 = v->getVectorBlock(0);  // V block
  Thyra::copy(*v0, x.ptr());

  return sStatus;
}

template <class Scalar>
void PhiEvaluator<Scalar>::applyMass(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> Mf,
                                     const Teuchos::RCP<const Thyra::VectorBase<Scalar>> f) const
{
  checkInitialized();
  phiLinSolv_->applyMass(Mf, f);
}

template <class Scalar>
void PhiEvaluator<Scalar>::solveMass(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> f,
                                     const Teuchos::RCP<const Thyra::VectorBase<Scalar>> Mf) const
{
  checkInitialized();
  phiLinSolv_->solveMass(f, Mf);
}

template <class Scalar>
void PhiEvaluator<Scalar>::applyJacobian(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> MJf,
                                         const Teuchos::RCP<const Thyra::VectorBase<Scalar>> f) const
{
  checkInitialized();
  phiLinSolv_->applyJacobian(MJf, f);
}


/*
 * PhiLinearSolver methods
 */
template <class Scalar>
void PhiLinearSolver<Scalar>::setLumpMassMatrix(bool lump)
{
  if (lumpMass_ != lump) {
    massMatrix_ == Teuchos::null;
    inverseMassMatrix_ == Teuchos::null;
  }

  lumpMass_ = lump;
}

template <class Scalar>
void PhiLinearSolver<Scalar>::setEigensolverParams(
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  eigensolverPL_ = pl;
}

template <class Scalar>
void PhiLinearSolver<Scalar>::initialize()
void PhiLinearSolver<Scalar>::clearMemory()
{
  massMatrix_ = Teuchos::null;
  inverseMassMatrix_ = Teuchos::null;
  jacobianMatrix_ = Teuchos::null;
}

template <class Scalar>
void PhiLinearSolver<Scalar>::checkInitialized(const PhiInitialization& mode) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      appModel_ == Teuchos::null, std::logic_error,
      "Error - PhiLinearSolver::initialize() Model not set!\n");
  TEUCHOS_TEST_FOR_EXCEPTION(
      massMatrix_ == Teuchos::null || inverseMassMatrix_ == Teuchos::null, std::logic_error,
      "Error - PhiLinearSolver::initialize() Mass matrix not computed!\n");

  if (mode == PhiInitialization::JACOBIAN_AND_MASS) {
    TEUCHOS_TEST_FOR_EXCEPTION(
        jacobianMatrix_ == Teuchos::null, std::logic_error,
        "Error - PhiLinearSolver::initialize() Jacobian matrix not computed!\n");
  }
}

template <class Scalar>
void PhiLinearSolver<Scalar>::computeMassMatrix(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs)
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef Thyra::ModelEvaluatorBase MEB;

  // first allocate space for the full mass matrix
  Teuchos::RCP<Thyra::LinearOpBase<Scalar>> fullMassMatrix = appModel_->create_W_op();

  // request only the mass matrix from the physics
  // Model evaluator builds: alpha*M*u_dot - beta*M*F(u) = 0,
  // where F(u) is the explicit tendency
  MEB::InArgs<Scalar> inArgs_new  = appModel_->createInArgs();
  inArgs_new.setArgs(inArgs);
  inArgs_new.set_x_dot(inArgs.get_x_dot());
  // Set x_dot to ensure we call the implicit model evaluator
  // for models that make a distiction based on x_dot==null
  // TODO: check this
  inArgs_new.set_alpha(1.0);
  inArgs_new.set_beta(0.0);

  //   TEUCHOS_TEST_FOR_EXCEPTION(inArgs_new.get_x().is_null(), std::runtime_error,
  //   "computeMassMatrix: x is null");

  // TEUCHOS_TEST_FOR_EXCEPTION(inArgs_new.get_x_dot().is_null(), std::runtime_error,
  //   "computeMassMatrix: x_dot is null");

  // set only the mass matrix
  MEB::OutArgs<Scalar> outArgs = appModel_->createOutArgs();
  outArgs.set_W_op(fullMassMatrix);

  // this will fill the mass matrix operator
  appModel_->evalModel(inArgs_new, outArgs);

  if(!lumpMass_) {
    massMatrix_ = fullMassMatrix;
    inverseMassMatrix_ = Thyra::inverse<Scalar>(*appModel_->get_W_factory(), fullMassMatrix);
  }
  else {
    // build lumped mass matrix (assumes all positive mass entries, does a simple sum)
    Teuchos::RCP<Thyra::VectorBase<Scalar>> ones = Thyra::createMember(*fullMassMatrix->domain());
    Thyra::assign(ones.ptr(), ST::one());

    Teuchos::RCP<Thyra::VectorBase<Scalar>> lumpedMassDiagonal = Thyra::createMember(*fullMassMatrix->range());
    Thyra::apply(*fullMassMatrix, Thyra::NOTRANS, *ones, lumpedMassDiagonal.ptr());

    //Teuchos::RCP<Thyra::VectorBase<Scalar> > invLumpedMassDiagonal = Thyra::createMember(*fullMassMatrix->range());
    Teuchos::RCP<Thyra::VectorBase<Scalar>> invLumpedMassDiagonal = ones;  // reuse memory from ones
    Thyra::reciprocal(*lumpedMassDiagonal, invLumpedMassDiagonal.ptr());

    massMatrix_ = Thyra::diagonal(lumpedMassDiagonal);
    inverseMassMatrix_ = Thyra::diagonal(invLumpedMassDiagonal);

    // if the mass matrix is lumped, keep the memory allocated for fullMassMatrix around for the Jacobian
    jacobianMatrix_ = fullMassMatrix;

    //Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    //lumpedMassDiagonal->describe(*out, Teuchos::VERB_EXTREME);
    //invLumpedMassDiagonal->describe(*out, Teuchos::VERB_EXTREME);
  }
  //Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  //massMatrix_->describe(*out, Teuchos::VERB_EXTREME);
  //inverseMassMatrix_->describe(*out, Teuchos::VERB_EXTREME);
}

template <class Scalar>
Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> PhiLinearSolver<Scalar>::buildL(const Scalar dt) const
{
  this->checkInitialized(PhiInitialization::JACOBIAN_AND_MASS);

  // Combine linear operators M_inv and J and multiply by -dt (minus is for implicit to explicit conversion)
  Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> L
    = Thyra::scale(-dt, Thyra::multiply<Scalar>(inverseMassMatrix_, jacobianMatrix_));

  return L;
}

template <class Scalar>
Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> PhiLinearSolver<Scalar>::buildATilde(
    const Scalar dt,
    const Teuchos::ArrayView<const Teuchos::RCP<const Thyra::VectorBase<Scalar>>> &rhs_B)
{
  this->checkInitialized(PhiInitialization::JACOBIAN_AND_MASS);

  const Thyra::Ordinal max_phi_order = rhs_B.size() - 1;
  auto KMatrix = this->buildK(max_phi_order);
  auto BMatrix = this->buildB(rhs_B);

  // Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  // invMassMatrix_->describe(*out, Teuchos::VERB_EXTREME);
  // jacobianMatrix_->describe(*out, Teuchos::VERB_EXTREME);

  // Combine linear operators M_inv and J and multiply by -dt (minus is for implicit to explicit conversion)
  Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> A
    = Thyra::scale(-dt, Thyra::multiply<Scalar>(inverseMassMatrix_, jacobianMatrix_));

  // Spaces
  auto V = A->domain();         // dim N, and also A->range()
  auto W = KMatrix->domain();   // dim max_phi_order, and also K->range()

  // Zero operator: V -> W  (for the (2,1) block)
  auto Z_VW = Thyra::zero<Scalar>(W, V);

  // Build block operator:
  // [ A  B ]
  // [ 0  K ]
  Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> Atilde = Thyra::block2x2<Scalar>(
      A,            // (1,1): V->V
      BMatrix,      // (1,2): W->V
      Z_VW,         // (2,1): V->W
      KMatrix       // (2,2): W->W
  );

  return Atilde;
}

template <class Scalar>
Teuchos::RCP<const Thyra::LinearOpBase<Scalar>>
PhiLinearSolver<Scalar>::buildK(const Thyra::Ordinal max_phi_order)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      max_phi_order < 1,
      std::invalid_argument,
      "buildK: max_phi_order must be positive.");

  // Space of dimension max_phi_order
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>> V =
      Thyra::defaultSpmdVectorSpace<Scalar>(max_phi_order);

  // Create a column multivector: rows = cols = max_phi_order
  Teuchos::RCP<Thyra::MultiVectorBase<Scalar>> K_mv = Thyra::createMembers(V, max_phi_order);

  // Initialize to zero
  Thyra::assign(K_mv.ptr(), Scalar(0.0));

  // Fill superdiagonal: K(i, i+1) = 1
  for (Thyra::Ordinal j = 1; j < max_phi_order; ++j)
  {
      // Column j, row j-1
      auto col_j = K_mv->col(j);
      Thyra::set_ele(j - 1, Scalar(1.0), col_j.ptr());
  }

  // Wrap as LinearOp and return
  return K_mv;
}

template <class Scalar>
Teuchos::RCP<const Thyra::LinearOpBase<Scalar>>
PhiLinearSolver<Scalar>::buildB(
  const Teuchos::ArrayView<const Teuchos::RCP<const Thyra::VectorBase<Scalar>>> &rhs_B)
{
  this->checkInitialized(PhiInitialization::ONLY_MASS);

  const Thyra::Ordinal max_phi_order = rhs_B.size() - 1;

  TEUCHOS_TEST_FOR_EXCEPTION(
      rhs_B.size() < 2,
      std::invalid_argument,
      "buildb: list of rhs must have at least two entries.");

  // N-dimensional space: use A's range (rows of an NxN operator)
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>> V_N = inverseMassMatrix_->range();

  // Create an N x max_phi_order multivector (N rows, max_phi_order columns)
  Teuchos::RCP<Thyra::MultiVectorBase<Scalar>> B_mv = Thyra::createMembers(V_N, max_phi_order);

  // Initialize to zero
  Thyra::assign(B_mv.ptr(), Scalar(0.0));

  // Fill the columns with rhs_B vectors in reverse order, excluding the first entry
  for (Thyra::Ordinal k = 0; k < max_phi_order; k++)
  {
    if (rhs_B[max_phi_order - k] != Teuchos::null)
    {
      auto col = B_mv->col(k);
      Thyra::assign(col.ptr(), *rhs_B[max_phi_order - k]);
    }
  }
  // wrap B_mv as LinOp and return
  return B_mv;
}

template <class Scalar>
Teuchos::RCP<Thyra::ProductVectorBase<Scalar>> PhiLinearSolver<Scalar>::buildv(
    const Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>> space,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>> x0) const
{
  // Create v and initialize to zero
  Teuchos::RCP<Thyra::ProductVectorBase<Scalar>> v =
    Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<Scalar>>(Thyra::createMember(space));

  // TODO: should the argument be a Thyra::ProductVectorSpace<Scalar>> instead of dynamic cast?

  Teuchos::RCP<Thyra::VectorBase<Scalar>> v0 = v->getNonconstVectorBlock(0);  // V block
  if (x0 != Teuchos::null)
  {
    Thyra::copy(*x0, v0.ptr());
  }
  else
  {
    Thyra::assign(v0.ptr(), Scalar(0.));
  }

  Teuchos::RCP<Thyra::VectorBase<Scalar>> v1 = v->getNonconstVectorBlock(1);  // W block
  Thyra::assign(v1.ptr(), Scalar(0.));

  // Get the last index
  const Thyra::Ordinal max_phi_order = v1->space()->dim();
  const Thyra::Ordinal g_last = max_phi_order - 1;

  // TODO: can this small dim x dim vector be distributed? It should live on every rank.
  // If this is a distributed space, set only on the owning rank
  if (auto spmdSpace = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<Scalar>>(space))
  {
    const Thyra::Ordinal localOffset = spmdSpace->localOffset();
    const Thyra::Ordinal localSubDim = spmdSpace->localSubDim();

    if (g_last >= localOffset && g_last < localOffset + localSubDim)
      Thyra::set_ele(g_last, Scalar(1), v1.ptr());
  }
  else
  {
    // Fallback for non-SPMD spaces
    Thyra::set_ele(g_last, Scalar(1), v1.ptr());
  }

  return v;
}

template <class Scalar>
void PhiLinearSolver<Scalar>::applyMass(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> Mf,
                                        const Teuchos::RCP<const Thyra::VectorBase<Scalar>> f) const
{
  this->checkInitialized(PhiInitialization::ONLY_MASS);
  // apply the mass matrix (either lumped, or not)
  if (f != Teuchos::null && Mf != Teuchos::null) {
    Thyra::apply(*massMatrix_, Thyra::NOTRANS, *f, Mf.ptr());
  }
}

template <class Scalar>
void PhiLinearSolver<Scalar>::solveMass(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> f,
                                        const Teuchos::RCP<const Thyra::VectorBase<Scalar>> Mf) const
{
  this->checkInitialized(PhiInitialization::ONLY_MASS);
  // invert the mass matrix (either lumped, or not)
  if (Mf != Teuchos::null && f != Teuchos::null) {
    Thyra::apply(*inverseMassMatrix_, Thyra::NOTRANS, *Mf, f.ptr());
  }
}

template <class Scalar>
void PhiLinearSolver<Scalar>::computeJacobian(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs)
{
  this->checkInitialized(PhiInitialization::ONLY_MASS);

  typedef Thyra::ModelEvaluatorBase MEB;

  // first allocate space for the Jacobian matrix
  jacobianMatrix_ = appModel_->create_W_op();

  // request only the Jacobian matrix from the physics
  // Model evaluator builds: alpha*u_dot - beta*M*F(u) = 0
  // where F(u) is the explicit tendency
  MEB::InArgs<Scalar> inArgs_new = appModel_->createInArgs();
  inArgs_new.setArgs(inArgs);
  inArgs_new.set_x_dot(inArgs.get_x_dot());
  // Set x_dot to ensure we call the implicit model evaluator
  // for models that make a distiction based on x_dot==null
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
  this->checkInitialized(PhiInitialization::JACOBIAN_AND_MASS);
  // apply the Jacobian matrix
  if (f != Teuchos::null && Jf != Teuchos::null) {
    Thyra::apply(*jacobianMatrix_, Thyra::NOTRANS, *f, Jf.ptr());
  }
}

template <class Scalar>
Thyra::SolveStatus<Scalar> PhiLinearSolver<Scalar>::assembleAndsolveMpJ(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
    const Teuchos::Ptr<Thyra::VectorBase<Scalar>> x,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>> Mf,
    Scalar alpha, Scalar beta) const
{
  // TODO: this method is not used and does not support mass lumping. Remove it.
  // However, it has some limited preconditioner support.
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  MEB::OutArgs<Scalar> outArgs = appModel_->createOutArgs();

  // first allocate space for the jacobian
  Teuchos::RCP<Thyra::LinearOpBase<Scalar>> MpJ = appModel_->create_W_op();
  Teuchos::RCP<Thyra::PreconditionerBase<Scalar>> MpJ_p = Teuchos::null;
  if (outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_W_prec)) {
    MpJ_p = appModel_->create_W_prec();
  }
  //TODO: support all types of preconditioner (compare NOX_Thyra_Group.C)

  MEB::InArgs<Scalar> inArgs_new = appModel_->createInArgs();
  inArgs_new.setArgs(inArgs);
  inArgs_new.set_x_dot(inArgs.get_x_dot());
  // Set x_dot to ensure we call the implicit model evaluator
  // for models that make a distiction based on x_dot==null
  // TODO: check this
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

  assign(x.ptr(), ST::zero());
  // Create solve criteria
  Thyra::SolveCriteria<Scalar> solveCriteria;
  solveCriteria.requestedTol = 1e-6; // p.get("Tolerance", 1.0e-6);

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

  // compute the solution to (MpJ)\Mf and write it to x
  Thyra::SolveStatus<Scalar> sStatus = LOWSB->solve(Thyra::NOTRANS, *Mf, x.ptr(), Teuchos::constPtr(solveCriteria));

  return sStatus;
}

template <class Scalar>
Thyra::SolveStatus<Scalar> PhiLinearSolver<Scalar>::solveMpJ(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> x,
                                                             const Teuchos::RCP<const Thyra::VectorBase<Scalar>> Mf,
                                                             Scalar alpha, Scalar beta) const
{
  this->checkInitialized(PhiInitialization::JACOBIAN_AND_MASS);

  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  // build an abstract MpJ linOp
  // TODO: figure out if it makes sense to add the physical matrices together, if both are the same Tpetra format.
  //       For lumpMass = true, this is probably fine as is.
  const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> aM = Thyra::scale<Scalar>(alpha, massMatrix_);
  const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> bJ = Thyra::scale<Scalar>(beta, jacobianMatrix_);
  const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> MpJ = Thyra::add(aM, bJ);

  // Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  // MpJ->describe(*out, Teuchos::VERB_EXTREME);

  // TODO: const-cast why?
  const Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar>> const_lowsFactory = appModel_->get_W_factory();
  const Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar>> lowsFactory =
    Teuchos::rcp_const_cast<Thyra::LinearOpWithSolveFactoryBase<Scalar>>(const_lowsFactory);

  // TODO: figure out how to add preconditioning
  const Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar>> LOWSB = Thyra::linearOpWithSolve(*lowsFactory, MpJ);

  assign(x.ptr(), ST::zero());
  // Create solve criteria
  Thyra::SolveCriteria<Scalar> solveCriteria;
  solveCriteria.requestedTol = 1e-6; //TODO: p.get("Tolerance", 1.0e-6);

  solveCriteria.solveMeasureType =
    Thyra::SolveMeasureType(Thyra::SOLVE_MEASURE_NORM_RESIDUAL, Thyra::SOLVE_MEASURE_NORM_INIT_RESIDUAL);

  // compute the solution to (MpJ)\Mf and write it to x
  Thyra::SolveStatus<Scalar> sStatus = LOWSB->solve(Thyra::NOTRANS, *Mf, x.ptr(), Teuchos::constPtr(solveCriteria));

  return sStatus;
}

template <class Scalar>
void PhiLinearSolver<Scalar>::computeJacobianSpectrumBounds(double& a, double& b, double& c)
{
  this->checkInitialized();

  using MT = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
  typedef Thyra::MultiVectorBase<Scalar> MV;
  typedef Thyra::LinearOpBase<Scalar>    OP;

  MT a_val =  Teuchos::ScalarTraits<MT>::rmax();
  MT b_val = -Teuchos::ScalarTraits<MT>::rmax();
  MT c_val =  Teuchos::ScalarTraits<MT>::zero();

  // Read eigensolver parameters from the stored ParameterList
  const int inev = eigensolverPL_->get<int>("Num Eigenvalues", 8);
  const int denseFallbackDim = eigensolverPL_->get<int>("Dense Fallback Dimension", 64);
  const int blockSize = eigensolverPL_->get<int>("Block Size", 2);
  const double eigTol = eigensolverPL_->get<double>("Convergence Tolerance", 0.08);
  const int maxRestarts = eigensolverPL_->get<int>("Maximum Restarts", 50);
  const bool relConvTol = eigensolverPL_->get<bool>("Relative Convergence Tolerance", true);
  const std::string which = eigensolverPL_->get<std::string>("Which", "LM");

  // Estimate up to the first nev eigenvalues
  const int dim = static_cast<int>(jacobianMatrix_->domain()->dim());
  const int nev      = std::min(inev, dim - 1);
  const int numBlocks = std::min(2 * nev, dim);

  // MinvJ = M^{-1}*J LinOp
  Teuchos::RCP<const OP> MinvJ =
      Thyra::scale<Scalar>(-1.0, Thyra::multiply<Scalar>(inverseMassMatrix_, jacobianMatrix_));

  // Fallback to dense eigen solve for small systems.
  if (dim <= denseFallbackDim) {
    // std::cout << "=== Fallback to Dense Eigen Solve ===" << std::endl;
    // Convert LinearOpBase to SerialDenseMatrix and compute eigenvalues directly.
    auto denseMinvJ = operatorToDense(MinvJ);
    Teuchos::Array<Scalar> eigs_re(dim);
    Teuchos::Array<Scalar> eigs_im(dim);
    denseEigenvalues(denseMinvJ, eigs_re, eigs_im);
    for (int i=0; i < eigs_re.size(); ++i) {
      auto evr = eigs_re[i];
      auto evi = eigs_im[i];
      // std::cout << "ev re: " << evr << " im: " << evi << std::endl;
      a_val = std::min(a_val, evr);
      b_val = std::max(b_val, evr);
      c_val = std::max(c_val, std::abs(evi));
    }
    a = std::min(0.0, a_val);
    b = std::max(0.0, b_val);
    c = c_val;
    return;
  }

#ifdef TEMPUS_ENABLE_ANASAZI

  TEUCHOS_TEST_FOR_EXCEPTION(
      (blockSize * numBlocks) <= nev, std::logic_error,
      "computeJacobianSpectrumBounds: requested number of eigenvalues (" << nev <<
      "). System dim (" << dim << "). The blockSize=" << blockSize << " or numBlocks=" << numBlocks <<
      " is too small. numBlocks*blockSize > nev is required.");

  // Init eigenvec multivector: warmstart from previous eigenvectors if available
  Teuchos::RCP<MV> initVec;
  const bool canWarmstart =
      evecsMinvJ_ != Teuchos::null &&
      static_cast<int>(evecsMinvJ_->domain()->dim()) == blockSize &&
      evecsMinvJ_->range()->isCompatible(*jacobianMatrix_->domain());

  if (canWarmstart) {
    // std::cout << "=== Warmstarting Krylov-Schur ===" << std::endl;
    initVec = evecsMinvJ_;
  } else {
    initVec = Thyra::createMembers(jacobianMatrix_->domain(), blockSize);
    Thyra::randomize(Scalar(-1), Scalar(1), initVec.ptr());
  }

  Teuchos::RCP<Anasazi::BasicEigenproblem<Scalar, MV, OP>> problem =
      Teuchos::rcp(new Anasazi::BasicEigenproblem<Scalar, MV, OP>());

  problem->setOperator(MinvJ);
  problem->setA(Thyra::scale<Scalar>(-1.0, jacobianMatrix_));
  problem->setM(massMatrix_);
  problem->setInitVec(initVec);
  problem->setNEV(nev);
  problem->setHermitian(false);  // M^{-1}*J is generally non-symmetric

  TEUCHOS_TEST_FOR_EXCEPTION(
      !problem->setProblem(), std::logic_error,
      "computeJacobianSpectrumBounds: BasicEigenproblem::setProblem() failed.");

  // Build the Anasazi solver parameter list from eigensolverPL_ with defaults.
  // "Num Blocks" is derived from nev and dim, so it is always computed here.
  // "Verbosity" is kept at Anasazi::Errors (not exposed in the user PL).
  Teuchos::ParameterList anasaziPL;
  anasaziPL.set("Which",                          which);
  anasaziPL.set("Block Size",                     blockSize);
  anasaziPL.set("Num Blocks",                     numBlocks);
  anasaziPL.set("Maximum Restarts",               maxRestarts);
  anasaziPL.set("Convergence Tolerance",          eigTol);
  anasaziPL.set("Relative Convergence Tolerance", relConvTol);
  anasaziPL.set("Verbosity",                      Anasazi::Errors);

  Anasazi::BlockKrylovSchurSolMgr<Scalar, MV, OP> solverMgr(problem, anasaziPL);
  solverMgr.solve();

  // Extract the current Ritz values. Can be unconverged
  const Anasazi::Eigensolution<Scalar, MV>& sol = problem->getSolution();
  std::vector<Anasazi::Value<Scalar>> evals;
  evals = solverMgr.getRitzValues();

  // Compute spectrum bounds
  // std::cout << "=== Krylov-Schur Ritz ===" << std::endl;
  for (const auto& ev : evals) {
    // std::cout << "ev re: " << ev.realpart << " im: " << ev.imagpart << std::endl;
    a_val = std::min(a_val, ev.realpart);
    b_val = std::max(b_val, ev.realpart);
    c_val = std::max(c_val, std::abs(ev.imagpart));
  }

  // Store converged eigenvectors for warmstart on the next call
  if (sol.numVecs > 0 && sol.Evecs != Teuchos::null) {
    const int cols = std::min(blockSize, sol.numVecs);
    evecsMinvJ_ = Thyra::createMembers(jacobianMatrix_->domain(), cols);
    Thyra::assign(evecsMinvJ_.ptr(),
                  *sol.Evecs->subView(Teuchos::Range1D(0, cols - 1)));
  }

  a = std::min(0.0, a_val);
  b = std::max(0.0, b_val);
  c = c_val;
#else
  TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error,
      "ANASAZI is unavailable and Block-Krylov-Schur Leja Ellipse update was requested. " <<
      "Reconfigure with Tempus_ENALBE_Anasazi=ON.\n");
#endif
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
