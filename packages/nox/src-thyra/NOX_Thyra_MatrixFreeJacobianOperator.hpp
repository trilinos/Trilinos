// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_THYRA_MATRIXFREE_JACOBIAN_OPERATOR_HPP
#define NOX_THYRA_MATRIXFREE_JACOBIAN_OPERATOR_HPP

#include "Thyra_LinearOpBase.hpp" // base class
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp" // base class

#include "NOX_Common.H"
#include "NOX_Utils.H"
#include "NOX_Solver_Generic.H"

// Forward Declarations

namespace NOX {
  namespace Abstract {
    class Group;
  }
}

namespace NOX {
namespace Thyra {

/*! \brief Concrete implementation of a Thyra::LinearOpBase object that approximates a Jacobian operator based on the Jacobian-Free Newton-Krylov method (see Knoll and Keyes Journal of Computational Physics 193 (2004) 357-397 for details).

This operator approximates the Jacobian action on a vector using rediual evaluations.  It is used by Newton-Krylov solvers since these methods only require the matrix-vector product \f$Jy\f$ in the linear solve iteration sequence.  This product can approximated by directional derivatives:

Forward (1st order accurate):
\f[ Jy \approx \frac{F(x + \delta y) - F(x)}{\delta} \f]

Central (2nd order accurate):
\f[ Jy \approx \frac{F(x + \delta y) - F(x - \delta y)}{\delta} \f]

where \f$J\f$ is the Jacobian, \f$y\f$ is the vector that the Jacobian is to be applied to, \f$x\f$ is the solution vector, \f$F(x)\f$ is the residual evaluation at \f$ x \f$,and \f$\delta\f$ is a scalar perturbation calculated by a variety of formulas (see MatrixFreeJacobianOperator::E_PerturbationType for references).

NOTE: The base evaluation \f$f(x)\f$ and its corresponding base solution \f$x\f$ is assumed to be supplied from an exterior source and that \f$f(x)\f$ is up to date wrt \f$x\f$.  Otherwise this algorithm would be inefficient due to extra base evaluations for every Jacobian-vector product.  A call must be made to either setBaseEvaluationToRawThyra() or setBaseEvaluationToNOXGroup() to set up the base objects.  If a NOX Group is used, the evaluation will be ensured to be correct because of the mechanics for isF().  In the case of raw objects, the user must take care to keep the base values up-to-date.
*/
template<typename Scalar>
class MatrixFreeJacobianOperator : public ::Thyra::LinearOpBase<Scalar>,
                   public Teuchos::ParameterListAcceptorDefaultBase {

 public:

  //! Define types for use of the perturbation parameter \f$ \delta\f$.
  enum E_DifferenceType {Forward=0, Centered=1};

  //! Defines the algorithm for computing .
  enum E_PerturbationType {
    //! Use the analogous algorithm by Salinger in LOCA v1.0 manual SAND2002-0396 p. 28 eqn. 2.43
    SalingerLOCA=0,
    //! Use algorithm developed for NOX by Kelley, Salinger and Pawlowski in 2001
    KelleySalingerPawlowski=1,
    //! Knoll and Keyes in JCP 193 (2004) 357-397
    KnollKeyes=2,
    //! Use a constant value defined by the user
    UserDefined=3
  };

  //! Defined where to get the base objects for the solution, \f$x\f$, and the residual, \f$ f(x) \f$, and the perturbed residual evaluation from.  In many cases the base values have been precomputed and can just be used to eliminate redundant residual evaluations.
  enum E_BaseEvaluationType
    {
      //! User defines thyra objects for solution vector, x, residual vector, f, and residual evaluations
      RawThyra,
      //! Use the nox group registered with this object
      NoxGroup
    };

  MatrixFreeJacobianOperator(Teuchos::ParameterList& printParams);

  /** \name Setup the base evaluation sources */
  //@{
  void setBaseEvaluationToRawThyra(const Teuchos::RCP<const ::Thyra::VectorBase<Scalar> >& x_base,
                   const Teuchos::RCP<const ::Thyra::VectorBase<Scalar> >& f_base,
                   const Teuchos::RCP< ::Thyra::ModelEvaluator<Scalar> > model);

  void setBaseEvaluationToNOXGroup(const Teuchos::RCP<const NOX::Abstract::Group>& base_group);
  //@}

  void setBaseInArgs(const Teuchos::RCP< ::Thyra::ModelEvaluatorBase::InArgs<Scalar> >& base_in_args);

  /** \name Derived from Teuchos::ParameterListAcceptor */
  //@{
  void setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& paramList);
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
  //@}

  /** \name Derived from Thyra::LinearOpBase */
  //@{

  Teuchos::RCP<const ::Thyra::VectorSpaceBase<Scalar> > range() const;

  Teuchos::RCP<const ::Thyra::VectorSpaceBase< Scalar > > domain() const;

  Teuchos::RCP<const ::Thyra::LinearOpBase< Scalar > > clone() const;

  bool opSupportedImpl (::Thyra::EOpTransp M_trans) const;

  void applyImpl(const ::Thyra::EOpTransp M_trans,
         const ::Thyra::MultiVectorBase< Scalar > &y,
         const Teuchos::Ptr< ::Thyra::MultiVectorBase< Scalar > > &u,
         const Scalar alpha,
         const Scalar beta) const;
  //@}

  //! Change the value of \f$ \lambda \f$ in the perturbation calculation.
  void setLambda(double lambda);

  //! Change the value of \f$ \lambda \f$ in the perturbation calculation.
  Scalar getLambda() const;

  //! Change the value of \f$ \delta \f$ in the perturbation calculation.
  void setUserDefinedDelta(double delta);

  //! Returns the user defined delta, \f$ \delta \f$, used for the perturbation.
  Scalar getUserDefinedDelta() const;

  //! Returns the last used value of delta \f$ \lambda \f$ in the perturbation calculation.
  Scalar getDelta() const;

protected:

  //! true if the algorithm has been setup using the setPArameterList() method
  bool setup_called_;

  //! User provided interface function
  Teuchos::RCP<const ::Thyra::ModelEvaluator<Scalar> > model_;

  //! Printing utilities.
  NOX::Utils utils_;

  //! A list of valid parameters
  mutable Teuchos::RCP<Teuchos::ParameterList> valid_params_;

  E_DifferenceType difference_type_;
  E_PerturbationType perturbation_type_;
  E_BaseEvaluationType base_evaluation_type_;

  //! Scale factor for eta calculation
  Scalar lambda_;

  //! Perturbation value to use in the directional derivative
  mutable Scalar delta_;

  //! Perturbation value to use in the directional derivative
  mutable Scalar user_defined_delta_;

  Teuchos::RCP<const ::Thyra::VectorBase<Scalar> > x_base_;

  Teuchos::RCP<const ::Thyra::VectorBase<Scalar> > f_base_;

  Teuchos::RCP< ::Thyra::VectorBase<Scalar> > x_perturb_;

  Teuchos::RCP< ::Thyra::VectorBase<Scalar> > f_perturb_;

  Teuchos::RCP< ::Thyra::VectorBase<Scalar> > f2_perturb_;

  Teuchos::RCP<const NOX::Thyra::Group> base_group_;

  mutable Teuchos::RCP< ::Thyra::ModelEvaluatorBase::InArgs<Scalar> > in_args_;

  mutable Teuchos::RCP< ::Thyra::ModelEvaluatorBase::OutArgs<Scalar> > out_args_;

};

} // namespace Thyra
} // namespace NOX

#include "NOX_Thyra_MatrixFreeJacobianOperator_impl.hpp"

#endif /* NOX_EPETRA_MATRIXFREE_JACOBIAN_OPERATOR_HPP */
