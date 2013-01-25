// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
       const std::vector< Teuchos::RCP<
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
      std::vector< Teuchos::RCP<
	LOCA::MultiContinuation::ConstraintInterfaceMVDX> > constraintMVDXPtrs;

      //! Composite constraint derivative
      Teuchos::RCP<NOX::Abstract::MultiVector> compositeDX;

    }; // Class CompositeConstraintMVDX

  } // namespace MultiContinuation

} // namespace LOCA

#endif // LOCA_MULTICONTINUATION_COMPOSITECONSTRAINT_H
