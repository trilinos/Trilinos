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

#ifndef LOCA_MULTICONTINUATION_COMPOSITECONSTRAINTMVDX_H
#define LOCA_MULTICONTINUATION_COMPOSITECONSTRAINTMVDX_H

#include "LOCA_MultiContinuation_ConstraintInterfaceMVDX.H" // base class
#include "LOCA_MultiContinuation_CompositeConstraint.H"     // base class

namespace LOCA { 

  namespace MultiContinuation { 

    /*! 
     * \brief Implementation of 
     * LOCA::MultiContinuation::ConstraintInterfaceMVDX for composite
     * constraints, i.e., a constraint comprised of multiple, separate
     * constraints.
     */
    class CompositeConstraintMVDX : 
      public virtual LOCA::MultiContinuation::ConstraintInterfaceMVDX,
      public virtual LOCA::MultiContinuation::CompositeConstraint {

    public:

      //! Constructor
      CompositeConstraintMVDX(
       const Teuchos::RCP<LOCA::GlobalData>& global_data,
       const vector< Teuchos::RCP<
       LOCA::MultiContinuation::ConstraintInterfaceMVDX> >& constraintObjects);
      
      //! Copy constructor
      CompositeConstraintMVDX(const CompositeConstraintMVDX& source, 
			      NOX::CopyType type = NOX::DeepCopy);

      //! Destructor
      ~CompositeConstraintMVDX();

      /*! 
       * @name Implementation of LOCA::MultiContinuation::ConstraintInterface 
       * virtual methods 
       */
      //@{

      //! Copy
      virtual void copy(const ConstraintInterface& source);

      //! Cloning function
      virtual 
      Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>
      clone(NOX::CopyType type = NOX::DeepCopy) const;

      //! Compute derivative of constraints w.r.t. solution vector x
      virtual NOX::Abstract::Group::ReturnType
      computeDX();

      //! Compute result_p = alpha * dg/dx * input_x
      /*!
       * Note that if there are n constraints and input_x has m columns,
       * result_p should be a n by m matrix and is equivalent to
       * \code
       *    input_x.multiply(alpha, dgdx, result_p).
       * \endcode
       */
      virtual NOX::Abstract::Group::ReturnType
      multiplyDX(double alpha, 
	      const NOX::Abstract::MultiVector& input_x,
	      NOX::Abstract::MultiVector::DenseMatrix& result_p) const;

      //! Compute result_x = alpha * dg/dx^T * op(b) + beta * result_x
      /*!
       * Note that this should be equivalent to
       * \code
       *    result_x.update(transb, alpha, dgdx, b, beta);
       * \endcode
       */
      virtual NOX::Abstract::Group::ReturnType
      addDX(Teuchos::ETransp transb,
	    double alpha, 
	    const NOX::Abstract::MultiVector::DenseMatrix& b,
	    double beta,
	    NOX::Abstract::MultiVector& result_x) const;

      //@}

      /*! 
       * @name Implementation of LOCA::MultiContinuation::ConstraintInterfaceMVDX 
       * virtual methods 
       */
      //@{

      //! Return solution component of constraint derivatives
      virtual const NOX::Abstract::MultiVector*
      getDX() const;

      //@}

    private:

      //! Prohibit generation and use of operator=()
      CompositeConstraintMVDX&
      operator=(const CompositeConstraintMVDX& source);

    protected:

      //! Array of constraintMVDX pointers
      vector< Teuchos::RCP<
	LOCA::MultiContinuation::ConstraintInterfaceMVDX> > constraintMVDXPtrs;

      //! Composite constraint derivative
      Teuchos::RCP<NOX::Abstract::MultiVector> compositeDX;

    }; // Class CompositeConstraintMVDX

  } // namespace MultiContinuation

} // namespace LOCA

#endif // LOCA_MULTICONTINUATION_COMPOSITECONSTRAINT_H
