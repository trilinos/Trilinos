// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
  template <class T> class ModelEvaluator;
  template <class T> class MultiVectorBase;
  template <class Scalar> class LinearOpWithSolveBase;
}

namespace NOX {

  namespace Thyra {

    //! A simple example of a group structure, based on Thyra.
    class Group : public virtual NOX::Abstract::Group {

    public:

      //! Constructor
      Group(const NOX::Thyra::Vector& initial_guess,
	    const Teuchos::RCP< const ::Thyra::ModelEvaluator<double> >& model);

      //! Copy constructor
      Group(const NOX::Thyra::Group& source, NOX::CopyType type = DeepCopy);

      //! Destructor.
      ~Group();

      NOX::Abstract::Group& operator=(const NOX::Abstract::Group& source);

      //! See above.
      NOX::Abstract::Group& operator=(const NOX::Thyra::Group& source);

      /** \brief . */
      Teuchos::RCP<const ::Thyra::VectorBase<double> > get_current_x() const;

      /** \brief . */
      Teuchos::RCP< ::Thyra::LinearOpWithSolveBase<double> >
      getNonconstJacobian();

      /** \brief . */
      Teuchos::RCP<const ::Thyra::LinearOpWithSolveBase<double> > getJacobian() const;

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

      //@}

      virtual Teuchos::RCP<NOX::Abstract::Group> 
      clone(NOX::CopyType type = NOX::DeepCopy) const;

      //! Print out the group
      void print() const;

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
      //@}
      
      //! Solution vector
      Teuchos::RCP<NOX::Thyra::Vector> x_vec_;
      
      //! Residual vector
      Teuchos::RCP<NOX::Thyra::Vector> f_vec_;
      
      //! Newton direction vector
      Teuchos::RCP<NOX::Thyra::Vector> newton_vec_;
      
      //! Gradient direction vector
      Teuchos::RCP<NOX::Thyra::Vector> gradient_vec_;
      
      //! Shared Jacobian operator
      Teuchos::RCP< 
        NOX::SharedObject<
          ::Thyra::LinearOpWithSolveBase<double>,
          NOX::Thyra::Group
          >
        > shared_jacobian_;
    
      //! Residual InArgs
      ::Thyra::ModelEvaluatorBase::InArgs<double> in_args_;
      
      //! Residual OutArgs
      ::Thyra::ModelEvaluatorBase::OutArgs<double> out_args_;

    };

  } // namespace LAPACK
} // namespace NOX


#endif
