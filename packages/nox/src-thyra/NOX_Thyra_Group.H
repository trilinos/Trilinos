// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef NOX_THYRA_GROUP_H
#define NOX_THYRA_GROUP_H

#include "Teuchos_RCP.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "NOX_Abstract_Group.H"	// base class
#include "NOX_Common.H"             // class data element (string)
#include "NOX_Thyra_Vector.H"	    // class data element
#include "NOX_SharedObjectTemplate.H"  // class data element 

// Forward declares
namespace NOX {
  namespace Parameter {
    class List;
  }
}

namespace Thyra {
  template <class Scalar> class ModelEvaluator;
  template <class Scalar> class MultiVectorBase;
  template <class Scalar> class LinearOpBase;
  template <class Scalar> class LinearOpWithSolveBase;
  template <class Scalar> class LinearOpWithSolveFactoryBase;
  template <class Scalar> class DefaultLinearOpSource;
  template <class Scalar> class PreconditionerFactoryBase;
  template <class Scalar> class PreconditionerBase;
}

namespace NOX {

  namespace Thyra {

    /** \brief A concrete implementation of the NOX::Abstract::Group using Thyra.

        NOTE: This Group supports row sum scaling of the function
        (residual and Jacobian).  This is enabled by setting a weight
        vector on the initial guess vector in the Group constructor.
        The residual and Jacobian must be scaled before and then
        unscaled after calls to construct the preconditioner and solve
        the linear system.  This follows the nox epetra group.  This
        design should be changed in a future nox refactor, but
        requires significant changes to the Group object.
    */
    class Group : public virtual NOX::Abstract::Group {

    public:

      /** \brief The default constructor that uses the linear solver from the ModelEvaluator.

          Most users should use this constructor.  It is meant to be
          used in conjunction with a stratimikos linear solver that is
          built as part of the input model evaluator.  For finer
          grained control over the use of the preconditioner and for
          Jacobian-Free Newton-Krylov cases, the power user
          constructor should be used.

          \param[in] initial_guess Initial guess for the solution vector
	  \param[in] model ModelEvaluator
	  \param[in] weight_vector Optional diagonal weighting vector for the model.
      */
      Group(const NOX::Thyra::Vector& initial_guess,
	    const Teuchos::RCP<const ::Thyra::ModelEvaluator<double> >& model,
	    const Teuchos::RCP<const ::Thyra::VectorBase<double> >& weight_vector = Teuchos::null);

      /** \brief Power user constructor that takes explicit linear solver objects to handle different combinations.

          This class allows the user to set user-defined linear
          operators and preconditioners (and corresponding
          factories). The user can set the linear_op to be a
          Jacobian-Free Newton Krylov operator (use the class
          NOX::Thyra::MatrixFreeJacobianOperator).

          \param[in] initial_guess Initial guess for the solution vector
	  \param[in] model ModelEvaluator
	  \param[in] linear_op Forward operator for the Jacobian.  Must be non-null
	  \param[in] lows_factory Factory for building and updating linear solver.
	  \param[in] prec_op Preconditioner operator.  If set to Teuchos::null, no preconditioner will be used.
	  \param[in] prec_factory Factory for updating the precondiitoner.  If set to Teuchos::null and there is a non-null prec_op, then the preconditioner will be updated using the model evaluator as long as the ModelEvaluator::outArgs supports W_prec.
	  \param[in] weight_vector Optional diagonal weighting vector for the model.
      */
      Group(const NOX::Thyra::Vector& initial_guess,
	    const Teuchos::RCP<const ::Thyra::ModelEvaluator<double> >& model,
	    const Teuchos::RCP< ::Thyra::LinearOpBase<double> >& linear_op,
	    const Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<double> >& lows_factory,
	    const Teuchos::RCP< ::Thyra::PreconditionerBase<double> >& prec_op,
	    const Teuchos::RCP< ::Thyra::PreconditionerFactoryBase<double> >& prec_factory,
	    const Teuchos::RCP<const ::Thyra::VectorBase<double> >& weight_vector = Teuchos::null);

      //! Copy constructor
      Group(const NOX::Thyra::Group& source, NOX::CopyType type = DeepCopy);

      //! Destructor.
      ~Group();

      NOX::Abstract::Group& operator=(const NOX::Abstract::Group& source);

      NOX::Abstract::Group& operator=(const NOX::Thyra::Group& source);

      Teuchos::RCP<const ::Thyra::VectorBase<double> > get_current_x() const;

      Teuchos::RCP< ::Thyra::LinearOpBase<double> >
      getNonconstJacobianOperator();

      Teuchos::RCP<const ::Thyra::LinearOpBase<double> >
      getJacobianOperator() const;

      Teuchos::RCP< ::Thyra::LinearOpWithSolveBase<double> >
      getNonconstJacobian();

      Teuchos::RCP<const ::Thyra::LinearOpWithSolveBase<double> > 
      getJacobian() const;

      /** @name "Compute" functions. */
      //@{

      void setX(const NOX::Abstract::Vector& y);
      //! See above
      void setX(const NOX::Thyra::Vector& y);

      void computeX(const NOX::Abstract::Group& grp, 
		    const NOX::Abstract::Vector& d, 
		    double step);
      //! See above.
      void computeX(const NOX::Thyra::Group& grp, 
		    const NOX::Thyra::Vector& d, 
		    double step);

      NOX::Abstract::Group::ReturnType computeF();

      NOX::Abstract::Group::ReturnType computeJacobian();

      NOX::Abstract::Group::ReturnType computeGradient();

      NOX::Abstract::Group::ReturnType 
      computeNewton(Teuchos::ParameterList& params);

      //@}

      /** @name Jacobian operations.
       *
       * Operations using the Jacobian matrix. These may not be defined in
       * matrix-free scenarios. */

      //@{
  
      NOX::Abstract::Group::ReturnType 
      applyJacobian(const NOX::Abstract::Vector& input, 
		    NOX::Abstract::Vector& result) const;
      
      NOX::Abstract::Group::ReturnType 
      applyJacobian(const NOX::Thyra::Vector& input, 
		    NOX::Thyra::Vector& result) const;

      NOX::Abstract::Group::ReturnType 
      applyJacobianMultiVector(const NOX::Abstract::MultiVector& input, 
			       NOX::Abstract::MultiVector& result) const;
      
      NOX::Abstract::Group::ReturnType 
      applyJacobianTranspose(const NOX::Abstract::Vector& input, 
			     NOX::Abstract::Vector& result) const;
      
      NOX::Abstract::Group::ReturnType 
      applyJacobianTranspose(const NOX::Thyra::Vector& input, 
			     NOX::Thyra::Vector& result) const;

      NOX::Abstract::Group::ReturnType 
      applyJacobianTransposeMultiVector(
				     const NOX::Abstract::MultiVector& input, 
				     NOX::Abstract::MultiVector& result) const;
      
      NOX::Abstract::Group::ReturnType 
      applyJacobianInverse(Teuchos::ParameterList& params, 
			   const NOX::Abstract::Vector& input, 
			   NOX::Abstract::Vector& result) const;

      NOX::Abstract::Group::ReturnType 
      applyJacobianInverse(Teuchos::ParameterList& params, 
			   const NOX::Thyra::Vector& input, 
			   NOX::Thyra::Vector& result) const;
      
      NOX::Abstract::Group::ReturnType 
      applyJacobianInverseMultiVector(
				    Teuchos::ParameterList& params, 
				    const NOX::Abstract::MultiVector& input, 
				    NOX::Abstract::MultiVector& result) const;
  
      //@}

      /** @name "Is" functions
       *
       * Checks to see if various objects have been computed. Returns true
       * if the corresponding "compute" function has been called since the
       * last update to the solution vector (via instantiation or
       * computeX). */
      //@{

      bool isF() const;
      bool isJacobian() const;
      bool isGradient() const;
      bool isNewton() const;

      //@}

      /** @name "Get" functions 
       *
       * Note that these function do not check whether or not the vectors
       * are valid. Must use the "Is" functions for that purpose. */
      //@{

      const NOX::Abstract::Vector& getX() const;

      const NOX::Abstract::Vector& getF() const;
  
      double getNormF() const;

      const NOX::Abstract::Vector& getGradient() const;

      const NOX::Abstract::Vector& getNewton() const;

      Teuchos::RCP< const NOX::Abstract::Vector > getXPtr() const;

      Teuchos::RCP< const NOX::Abstract::Vector > getFPtr() const;
  
      Teuchos::RCP< const NOX::Abstract::Vector > getGradientPtr() const;

      Teuchos::RCP< const NOX::Abstract::Vector > getNewtonPtr() const;

      //@}

      virtual Teuchos::RCP<NOX::Abstract::Group> 
      clone(NOX::CopyType type = NOX::DeepCopy) const;

      //! Print out the group
      void print() const;

      /** FOR POWER USERS ONLY!  Grab the inargs used by nox and allow the user to change it.  Used by pseudo-transient solver to add an x_dot and alpha/beta to the model evaluator call. */
      ::Thyra::ModelEvaluatorBase::InArgs<double>& getNonconstInArgs();

      const ::Thyra::ModelEvaluatorBase::InArgs<double>& getInArgs() const;

      Teuchos::RCP< const ::Thyra::ModelEvaluator<double> > getModel() const;

    protected:

      //! resets the isValid flags to false
      void resetIsValidFlags();

      //! Apply Jacobian inverse using Thyra objects
      NOX::Abstract::Group::ReturnType 
      applyJacobianInverseMultiVector(
			      Teuchos::ParameterList& p, 
			      const ::Thyra::MultiVectorBase<double>& input, 
			      ::Thyra::MultiVectorBase<double>& result) const;

      ::Thyra::ESolveMeasureNormType 
      getThyraNormType(const string& name) const;

      //! Finalizes LOWS to be a valid solver for the Jacobian
      void updateLOWS() const;

      void scaleResidualAndJacobian() const;

      void unscaleResidualAndJacobian() const;

    protected:

      //! Problem interface
      Teuchos::RCP< const ::Thyra::ModelEvaluator<double> > model_;
      
      /** @name IsValid flags 
       *  
       * True if the current solution is up-to-date with respect to the
       * currect solution vector. */
      //@{
      bool is_valid_f_;
      bool is_valid_jacobian_;
      bool is_valid_newton_dir_;
      bool is_valid_gradient_dir_;
      mutable bool is_valid_lows_;
      //@}
      
      //! Solution vector
      Teuchos::RCP<NOX::Thyra::Vector> x_vec_;
      
      //! Residual vector
      Teuchos::RCP<NOX::Thyra::Vector> f_vec_;
      
      //! Newton direction vector
      Teuchos::RCP<NOX::Thyra::Vector> newton_vec_;
      
      //! Gradient direction vector
      Teuchos::RCP<NOX::Thyra::Vector> gradient_vec_;
      
      //! Shared Jacobian operator with solve
      Teuchos::RCP< 
        NOX::SharedObject<
          ::Thyra::LinearOpWithSolveBase<double>,
          NOX::Thyra::Group
          >
        > shared_jacobian_;
    
      //! Jacobian operator 
      Teuchos::RCP< ::Thyra::LinearOpBase<double> > lop_;

      //! Thyra LOWS factory for building Jacobians
      Teuchos::RCP< const ::Thyra::LinearOpWithSolveFactoryBase<double> > lows_factory_;

      //! Source base needed to create preconditioner
      Teuchos::RCP< const ::Thyra::DefaultLinearOpSource<double> > losb_;

      //! Preconditioner for Jacobian
      Teuchos::RCP< ::Thyra::PreconditionerBase<double> > prec_;

      //! Preconditioner factory
      Teuchos::RCP< ::Thyra::PreconditionerFactoryBase<double> > prec_factory_;

      //! Residual InArgs
      mutable ::Thyra::ModelEvaluatorBase::InArgs<double> in_args_;
      
      //! Residual OutArgs
      mutable ::Thyra::ModelEvaluatorBase::OutArgs<double> out_args_;

      /** \brief Optional wieghting vector for function scaling

           This is pulled out of the initial guess vector
      */
      Teuchos::RCP<const ::Thyra::VectorBase<double> > weight_vec_;

      //! Inverse of weight vector used to unscale problem
      mutable Teuchos::RCP< ::Thyra::VectorBase<double> > inv_weight_vec_;

    };

  } // namespace LAPACK
} // namespace NOX


#endif
