// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_BELOS_LINEAR_OP_WITH_SOLVE_FACTORY_DECL_HPP
#define THYRA_BELOS_LINEAR_OP_WITH_SOLVE_FACTORY_DECL_HPP

#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"


namespace Teuchos { class ParameterList; }


namespace Thyra {


//! \brief
enum EBelosSolverType {
  SOLVER_TYPE_BLOCK_GMRES,
  SOLVER_TYPE_PSEUDO_BLOCK_GMRES,
  SOLVER_TYPE_BLOCK_CG,
  SOLVER_TYPE_PSEUDO_BLOCK_CG,
  SOLVER_TYPE_PSEUDO_BLOCK_STOCHASTIC_CG,
  SOLVER_TYPE_GCRODR,
  SOLVER_TYPE_RCG,
  SOLVER_TYPE_MINRES,
  SOLVER_TYPE_TFQMR,
  SOLVER_TYPE_BICGSTAB,
  SOLVER_TYPE_FIXEDPOINT,
  SOLVER_TYPE_TPETRA_GMRES,
  SOLVER_TYPE_TPETRA_GMRES_PIPELINE,
  SOLVER_TYPE_TPETRA_GMRES_SINGLE_REDUCE,
  SOLVER_TYPE_TPETRA_GMRES_SSTEP
};


//! \brief
inline std::istream& operator>>(
  std::istream& is, EBelosSolverType& sType)
{
  int intval;
  is >> intval;
  sType = (EBelosSolverType)intval;
  return is;
}


/** \brief <tt>LinearOpWithSolveFactoryBase</tt> subclass implemented in terms
 * of <tt>Belos</tt>.
 *
 * ToDo: Finish Documentation!
 *
 * \ingroup Belos_Thyra_adapters_grp
 */
template<class Scalar>
class BelosLinearOpWithSolveFactory : public LinearOpWithSolveFactoryBase<Scalar> {
public:

  /** \name Public types */
  //@{
  /** \brief . */

  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType  MagnitudeType;

  //@}

  /** \name Parameter names for Parameter List */
  //@{

  /** \brief . */
  static const std::string  SolverType_name;
  /** \brief . */           
  static const std::string  SolverType_default;
  /** \brief . */
  static const std::string  SolverTypes_name;
  /** \brief . */
  static const std::string  BlockGMRES_name;
  /** \brief . */
  static const std::string  PseudoBlockGMRES_name;
  /** \brief . */
  static const std::string  BlockCG_name;
  /** \brief . */
  static const std::string  PseudoBlockCG_name;
  /** \brief . */
  static const std::string  PseudoBlockStochasticCG_name;
  /** \brief . */
  static const std::string  GCRODR_name;
  /** \brief . */
  static const std::string  RCG_name;
  /** \brief . */
  static const std::string  MINRES_name;
  /** \brief . */
  static const std::string  TFQMR_name;
  /** \brief . */
  static const std::string  BiCGStab_name;
  /** \brief . */
  static const std::string  FixedPoint_name;
  /** \brief . */
  static const std::string  TpetraGmres_name;
  /** \brief . */
  static const std::string  TpetraGmresPipeline_name;
  /** \brief . */
  static const std::string  TpetraGmresSingleReduce_name;
  /** \brief . */
  static const std::string  TpetraGmresSstep_name;
  /** \brief . */
  static const std::string  ConvergenceTestFrequency_name;

  //@}

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Construct without preconditioner factory. */
  BelosLinearOpWithSolveFactory();

  /** \brief Calls <tt>this->setPreconditionerFactory(precFactory)</tt.  . */
  BelosLinearOpWithSolveFactory(
    const Teuchos::RCP<PreconditionerFactoryBase<Scalar> >      &precFactory
    );

  //@}

  /** @name Overridden public functions from LinearOpWithSolveFactoryBase */
  //@{
  /** \brief Returns true . */
  bool acceptsPreconditionerFactory() const;
  /** \brief . */
  void setPreconditionerFactory(
    const Teuchos::RCP<PreconditionerFactoryBase<Scalar> >      &precFactory
    ,const std::string                                          &precFactoryName
    );
  /** \brief . */
  Teuchos::RCP<PreconditionerFactoryBase<Scalar> > getPreconditionerFactory() const;
  /** \brief . */
  void unsetPreconditionerFactory(
    Teuchos::RCP<PreconditionerFactoryBase<Scalar> >            *precFactory
    ,std::string                                                *precFactoryName
    );
  /** \brief . */
  bool isCompatible( const LinearOpSourceBase<Scalar> &fwdOpSrc ) const;
  /** \brief . */
  Teuchos::RCP<LinearOpWithSolveBase<Scalar> > createOp() const;
  /** \brief . */
  void initializeOp(
    const Teuchos::RCP<const LinearOpSourceBase<Scalar> >       &fwdOpSrc
    ,LinearOpWithSolveBase<Scalar>                              *Op
    ,const ESupportSolveUse                                     supportSolveUse
    ) const;
  /** \brief . */
  void initializeAndReuseOp(
    const Teuchos::RCP<const LinearOpSourceBase<Scalar> >       &fwdOpSrc
    ,LinearOpWithSolveBase<Scalar>                              *Op
    ) const;
  /** \brief . */
  void uninitializeOp(
    LinearOpWithSolveBase<Scalar>                               *Op
    ,Teuchos::RCP<const LinearOpSourceBase<Scalar> >            *fwdOpSrc
    ,Teuchos::RCP<const PreconditionerBase<Scalar> >            *prec
    ,Teuchos::RCP<const LinearOpSourceBase<Scalar> >            *approxFwdOpSrc
    ,ESupportSolveUse                                           *supportSolveUse
    ) const;
  /** \brief . */
  bool supportsPreconditionerInputType(const EPreconditionerInputType precOpType) const;
  /** \brief . */
  void initializePreconditionedOp(
    const Teuchos::RCP<const LinearOpSourceBase<Scalar> >       &fwdOpSrc
    ,const Teuchos::RCP<const PreconditionerBase<Scalar> >      &prec
    ,LinearOpWithSolveBase<Scalar>                              *Op
    ,const ESupportSolveUse                                     supportSolveUse
    ) const;
  /** \brief . */
  void initializeApproxPreconditionedOp(
    const Teuchos::RCP<const LinearOpSourceBase<Scalar> >       &fwdOpSrc
    ,const Teuchos::RCP<const LinearOpSourceBase<Scalar> >      &approxFwdOpSrc
    ,LinearOpWithSolveBase<Scalar>                              *Op
    ,const ESupportSolveUse                                     supportSolveUse
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
  // Private types


  // /////////////////////////
  // Private data members

  Teuchos::RCP<PreconditionerFactoryBase<Scalar> >  precFactory_;
  std::string                                       precFactoryName_;
  Teuchos::RCP<Teuchos::ParameterList>              thisValidParamList_;
  Teuchos::RCP<Teuchos::ParameterList>              paramList_;
  EBelosSolverType solverType_;
  int convergenceTestFrequency_;

  // /////////////////////////
  // Private member functions

  static Teuchos::RCP<const Teuchos::ParameterList> generateAndGetValidParameters();

  void updateThisValidParamList();

  void initializeOpImpl(
    const Teuchos::RCP<const LinearOpSourceBase<Scalar> >       &fwdOpSrc
    ,const Teuchos::RCP<const LinearOpSourceBase<Scalar> >      &approxFwdOpSrc
    ,const Teuchos::RCP<const PreconditionerBase<Scalar> >      &prec
    ,const bool                                                 reusePrec
    ,LinearOpWithSolveBase<Scalar>                              *Op
    ,const ESupportSolveUse                                     supportSolveUse
    ) const;

};

//@}

} // namespace Thyra

#endif // THYRA_BELOS_LINEAR_OP_WITH_SOLVE_FACTORY_DECL_HPP
