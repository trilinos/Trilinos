// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_AMESOS2_LINEAR_OP_WITH_SOLVE_FACTORY_DECL_HPP
#define THYRA_AMESOS2_LINEAR_OP_WITH_SOLVE_FACTORY_DECL_HPP

#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_Amesos2Types.hpp"
#include "Amesos2.hpp"
#include "Amesos2_Version.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Trilinos_Details_LinearSolver.hpp"
#include "Trilinos_Details_LinearSolverFactory.hpp"
#include "Thyra_Amesos2LinearOpWithSolve_decl.hpp"

namespace Thyra {

template<typename Scalar>
class Amesos2LinearOpWithSolveFactory : public LinearOpWithSolveFactoryBase<Scalar> {
public:
  using ConverterT = typename Amesos2LinearOpWithSolve<Scalar>::ConverterT;
  using MAT = typename Amesos2LinearOpWithSolve<Scalar>::MAT;
  using MV = typename Amesos2LinearOpWithSolve<Scalar>::MV;
  using Solver = typename Amesos2LinearOpWithSolve<Scalar>::Solver;

  /** \name Parameter names for Paramter List */
  //@{

  /** \brief . */
  static const std::string SolverType_name;
  /** \brief . */
  static const std::string RefactorizationPolicy_name;
  /** \brief . */
  static const std::string ThrowOnPreconditionerInput_name;
  /** \brief . */
  static const std::string Amesos2_Settings_name;

  //@}

  /** \name Constructors/initializers/accessors */
  //@{

  /** \brief . */
  ~Amesos2LinearOpWithSolveFactory();

  /** \brief Constructor which sets the defaults.
   */
  Amesos2LinearOpWithSolveFactory(
    const Amesos2::ESolverType solverType
#ifdef HAVE_AMESOS2_KLU2
    = Amesos2::KLU2,
#else
    = Amesos2::LAPACK,
#endif
    const Amesos2::ERefactorizationPolicy refactorizationPolicy
      = Amesos2::REPIVOT_ON_REFACTORIZATION,
    const bool throwOnPrecInput = true
    );
    
  //@}

  /** @name Overridden public functions from LinearOpWithSolveFactoryBase */
  //@{

  /** \brief . */
  bool isCompatible( const LinearOpSourceBase<Scalar> &fwdOpSrc ) const;

  /** \brief . */
  Teuchos::RCP<LinearOpWithSolveBase<Scalar> > createOp() const;

  /** \brief . */
  void initializeOp(
    const Teuchos::RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
    LinearOpWithSolveBase<Scalar> *Op,
    const ESupportSolveUse supportSolveUse
    ) const;

  /** \brief Returns <tt>false</tt> . */
  bool supportsPreconditionerInputType(const EPreconditionerInputType precOpType) const;

  /** \brief Throws exception if <tt>this->throwOnPrecInput()==true</tt> and
   * calls <tt>this->initializeOp(fwdOpSrc,Op)</tt> otherwise
   */
  void initializePreconditionedOp(
    const Teuchos::RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
    const Teuchos::RCP<const PreconditionerBase<Scalar> > &prec,
    LinearOpWithSolveBase<Scalar> *Op,
    const ESupportSolveUse supportSolveUse
    ) const;

  /** \brief Throws exception if <tt>this->throwOnPrecInput()==true</tt> and
   * calls <tt>this->initializeOp(fwdOpSrc,Op)</tt> otherwise
   */
  void initializePreconditionedOp(
    const Teuchos::RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
    const Teuchos::RCP<const LinearOpSourceBase<Scalar> > &approxFwdOpSrc,
    LinearOpWithSolveBase<Scalar> *Op,
    const ESupportSolveUse supportSolveUse
    ) const;

  /** \brief . */
  void uninitializeOp(
    LinearOpWithSolveBase<Scalar> *Op,
    Teuchos::RCP<const LinearOpSourceBase<Scalar> > *fwdOpSrc,
    Teuchos::RCP<const PreconditionerBase<Scalar> > *prec,
    Teuchos::RCP<const LinearOpSourceBase<Scalar> > *approxFwdOpSrc,
    ESupportSolveUse *supportSolveUse
    ) const;

  //@}

  /** @name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList);
  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const;
  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  //@}

  /** \name Public functions overridden from Teuchos::Describable. */
  //@{

  /** \brief . */
  std::string description() const;

  //@}

private:

  // /////////////////////////
  // Private data members

  Amesos2::ESolverType solverType_;
  Amesos2::ERefactorizationPolicy refactorizationPolicy_;
  bool throwOnPrecInput_;
  Teuchos::RCP<Teuchos::ParameterList> paramList_;

  // /////////////////////////
  // Private member functions

  static Teuchos::RCP<const Teuchos::ParameterList>
  generateAndGetValidParameters();

};

} // namespace Thyra

#endif // THYRA_AMESOS2_LINEAR_OP_WITH_SOLVE_FACTORY_HPP
