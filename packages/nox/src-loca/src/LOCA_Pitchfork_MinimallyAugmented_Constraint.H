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

#ifndef LOCA_PITCHFORK_MINIMALLYAUGMENTED_CONSTRAINT_H
#define LOCA_PITCHFORK_MINIMALLYAUGMENTED_CONSTRAINT_H

#include "LOCA_TurningPoint_MinimallyAugmented_Constraint.H" // base class

// forward declarations
namespace LOCA {
  namespace Pitchfork {
    namespace MinimallyAugmented {
      class AbstractGroup;
    }
  }
}

namespace LOCA { 

  namespace Pitchfork {

    namespace MinimallyAugmented {

      /*! 
       * \brief Implementation of 
       * LOCA::MultiContinuation::ConstraintInterfaceMVDX for computing
       * pitchforks for the minimally augmented pitchfork formulation.
       */
      /*!
       * This class implements the pitchfork constraint equations
       * \f$\sigma(x,p) = 0\f$, \f$\langle \psi,x \rangle = 0\f$ for 
       * the minimally augmented pitchfork formulation where \f$\sigma\f$ 
       * is defined via
       * \f[
       *     \begin{bmatrix}
       *       J   & a \\
       *       b^T & 0 
       *     \end{bmatrix}
       *     \begin{bmatrix}
       *       v \\
       *       \sigma_1
       *     \end{bmatrix} = 
       *     \begin{bmatrix}
       *       0 \\
       *       n
       *     \end{bmatrix},
       * \f]
       * \f[
       *     \begin{bmatrix}
       *       J^T & b \\
       *       a^T & 0 
       *     \end{bmatrix}
       *     \begin{bmatrix}
       *       w \\
       *       \sigma_2
       *     \end{bmatrix} = 
       *     \begin{bmatrix}
       *       0 \\
       *       n
       *     \end{bmatrix},
       * \f]
       * \f[
       *     \sigma = -w^T J v/n
       * \f]
       * for any vectors \f$a\f$ and \f$b\f$ in \f$\Re^n\f$.  Using these
       * relationships, it is easy to show
       * \f[
       *     \begin{split}
       *        \sigma_x &= -(w^T J v)_x/n = -w^T J_x v/n \\
       *        \sigma_p &= -(w^T J v)_p/n = -w^T J_p v/n
       *     \end{split}
       * \f]
       *
       * The class is derived from 
       * LOCA::TurningPoint::MinimallyAugmented::Constraint.  See this class 
       * for a description of available parameters.
       */
      class Constraint : 
	public LOCA::TurningPoint::MinimallyAugmented::Constraint {

      public:

	//! Constructor
	Constraint(
	 const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	 const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<Teuchos::ParameterList>& pfParams,
	 const Teuchos::RefCountPtr<LOCA::Pitchfork::MinimallyAugmented::AbstractGroup>& g,
	 bool is_symmetric,
	 const NOX::Abstract::Vector& a,
	 const NOX::Abstract::Vector* b,
	 const Teuchos::RefCountPtr<const NOX::Abstract::Vector>& psi,
	 int bif_param);

	//! Copy constructor
	Constraint(const Constraint& source, 
		   NOX::CopyType type = NOX::DeepCopy);

	//! Destructor
	virtual ~Constraint();

	//! Set the group pointer
	/*!
	 * This method should be called when ever the constrained group
	 * is copied, since we don't explicitly copy the underlying group
	 * here.
	 */
	virtual void setGroup(const Teuchos::RefCountPtr<LOCA::TurningPoint::MinimallyAugmented::AbstractGroup>& g);

	/*! 
	 * @name Implementation of LOCA::MultiContinuation::ConstraintInterface
	 * virtual methods 
	 */
	//@{

	//! Copy
	virtual void 
	copy(const LOCA::MultiContinuation::ConstraintInterface& source);

	//! Cloning function
	virtual 
	Teuchos::RefCountPtr<LOCA::MultiContinuation::ConstraintInterface>
	clone(NOX::CopyType type = NOX::DeepCopy) const;

	//! Return number of constraints
	virtual int numConstraints() const;

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

	//! Return constraint residuals
	virtual const NOX::Abstract::MultiVector::DenseMatrix&
	getConstraints() const;

	//! Return solution component of constraint derivatives
	virtual const NOX::Abstract::MultiVector*
	getDX() const;

	//@}

      private:

	//! Prohibit generation and use of operator=()
	Constraint& operator=(const Constraint& source);

      protected:

	//! Pointer to pitchfork group
	Teuchos::RefCountPtr<LOCA::Pitchfork::MinimallyAugmented::AbstractGroup> pf_grp;

	//! Vector for \f$\psi\f$
	Teuchos::RefCountPtr<const NOX::Abstract::Vector> psi_vector;

	//! Stores dg/dx
	Teuchos::RefCountPtr<NOX::Abstract::MultiVector> dgdx;

	//! Constraint values
	NOX::Abstract::MultiVector::DenseMatrix pf_constraints;

      }; // Class Constraint

    } // namespace MinimallyAugmented

  } // namespace Pitchfork

} // namespace LOCA

#endif // LOCA_PITCHFORK_MINIMALLYAUGMENTED_CONSTRAINT_H
