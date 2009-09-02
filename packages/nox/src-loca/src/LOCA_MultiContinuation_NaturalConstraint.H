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

#ifndef LOCA_MULTICONTINUATION_NATURALCONSTRAINT_H
#define LOCA_MULTICONTINUATION_NATURALCONSTRAINT_H

#include "LOCA_MultiContinuation_ConstraintInterfaceMVDX.H" // base class
#include "NOX_Abstract_MultiVector.H"                   // for dense matrix

// forward declarations
namespace LOCA {
  class GlobalData;
  namespace MultiContinuation {
    class NaturalGroup;
  }
}

namespace LOCA { 

  namespace MultiContinuation { 

    /*! 
     * \brief Implementation of LOCA::MultiContinuation::ConstraintInterface
     * for natural continuation.
     */
    /*!
     * This class implements the natural constraint equation for
     * natural continuation:
     * \f[
     *     g(x,p,x_0,p_0,x^\ast,p^\ast,v,\Delta s)=
     *         p-p_0-v_p \Delta s
     * \f]
     * where \f$v_p\f$ is the parameter component of the predictor direction 
     * \f$v\f$.
     */
    class NaturalConstraint : 
      public LOCA::MultiContinuation::ConstraintInterfaceMVDX {

    public:

      //! Constructor
      NaturalConstraint(
     const Teuchos::RCP<LOCA::GlobalData>& global_data,
     const Teuchos::RCP<LOCA::MultiContinuation::NaturalGroup>& grp);

      //! Copy constructor
      NaturalConstraint(const NaturalConstraint& source, 
			NOX::CopyType type = NOX::DeepCopy);

      //! Destructor
      ~NaturalConstraint();

      //! Set pointer to natural group
      virtual void setNaturalGroup(const Teuchos::RCP<LOCA::MultiContinuation::NaturalGroup>& grp);

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

      //! Return number of constraints
      virtual int numConstraints() const;

      //! Set the solution vector to y.
      virtual void setX(const NOX::Abstract::Vector& y);

      //! Sets parameter indexed by paramID
      virtual void setParam(int paramID, double val);

      //! Sets parameters indexed by paramIDs
      virtual void setParams(
			 const vector<int>& paramIDs, 
			 const NOX::Abstract::MultiVector::DenseMatrix& vals);

      //! Compute continuation constraint equations
      virtual NOX::Abstract::Group::ReturnType
      computeConstraints();

      //! Compute derivative of constraints w.r.t. solution vector x
      virtual NOX::Abstract::Group::ReturnType
      computeDX();

      //! Compute derivative of constraints w.r.t. supplied parameters.
      /*!
       * The first column of \c dgdp should be filled with the constraint
       * residuals \f$g\f$ if \c isValidG is \c false.  If \c isValidG is
       * \c true, then the \c dgdp contains \f$g\f$ on input.
       */
      virtual NOX::Abstract::Group::ReturnType
      computeDP(const vector<int>& paramIDs, 
		NOX::Abstract::MultiVector::DenseMatrix& dgdp, 
		bool isValidG);

      //! Return \c true if constraint residuals are valid
      virtual bool isConstraints() const;

      //! Return \c true if derivatives of constraints w.r.t. x are valid
      virtual bool isDX() const;

      //! Return constraint residuals
      virtual const NOX::Abstract::MultiVector::DenseMatrix&
      getConstraints() const;

      //! Return solution component of constraint derivatives
      /*!
       * Since the solution component of the derivative is always zero,
       * this always returns NULL.
       */
      virtual const NOX::Abstract::MultiVector*
      getDX() const;

      /*! 
       * \brief Return \c true if solution component of constraint 
       * derivatives is zero
       */
      virtual bool isDXZero() const;

      //@}

    private:

      //! Prohibit generation and use of operator=()
      NaturalConstraint& operator=(const NaturalConstraint& source);

    protected:

      //! LOCA global data object
      Teuchos::RCP<LOCA::GlobalData> globalData;

      //! Pointer to arc-length group
      Teuchos::RCP<
	LOCA::MultiContinuation::NaturalGroup> naturalGroup;

      //! Constraint values
      NOX::Abstract::MultiVector::DenseMatrix constraints;

      //! Flag indicating whether constraints are valid
      bool isValidConstraints;

      //! Continuation parameter IDs
      vector<int> conParamIDs;

    }; // Class ArcLengthConstraint

  } // namespace MultiContinuation

} // namespace LOCA

#endif // LOCA_MULTICONTINUATION_ARCLENGTHCONSTRAINT_H
