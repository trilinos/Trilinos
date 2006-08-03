// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//                  LOCA Continuation Algorithm Package
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov) or Eric Phipps
// (etphipp@sandia.gov), Sandia National Laboratories.
//
// ************************************************************************
//@HEADER

#ifndef LOCA_TURNINGPOINT_MOORESPENCE_SALINGERBORDERING_H
#define LOCA_TURNINGPOINT_MOORESPENCE_SALINGERBORDERING_H

#include "LOCA_TurningPoint_MooreSpence_SolverStrategy.H"
#include "NOX_Abstract_MultiVector.H"

// forward declarations
namespace LOCA {
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }
}

namespace LOCA {

  namespace TurningPoint {

    namespace MooreSpence {
 
      /*!
       * \brief Moore-Spence turning point solver strategy based on
       * "Salinger" bordering.  This is the classic 4-solve bordering
       * method.
       */
      /*!
       * This class solves the Moore-Spence turning point Newton equations:
       * \f[
       *     \begin{bmatrix}
       *          J   &    0   &   f_p \\
       *       (Jv)_x &    J   & (Jv)_p \\
       *          0   & \phi^T &   0
       *     \end{bmatrix}
       *     \begin{bmatrix}
       *        X \\
       *        Y \\
       *        z
       *     \end{bmatrix} = 
       *     \begin{bmatrix}
       *        F \\
       *        G \\
       *        h
       *     \end{bmatrix}
       *  \f]
       *  via the following block elimination scheme:
       *  \f[
       *    \begin{split}
       *       J [A \; b] &= [F \; f_p] \\
       *       J [C \; d] &= (Jv)_x[A \; b] - [G \; (Jv)_p] \\
       *       z &= (h + \phi^T C) / \phi^T d \\
       *       X &= A - b z \\
       *       Y &= -C + d z
       *    \end{split}
       *   \f]
       */
      class SalingerBordering : 
	public LOCA::TurningPoint::MooreSpence::SolverStrategy {

      public:

	//! Constructor. 
	/*!
	 * \param global_data [in] Global data object
	 * \param topParams [in] Parsed top-level parameter list 
	 * \param solverParams [in] Bordered solver parameters.  Currently none
	 * are referenced.
	 */
	SalingerBordering(
	 const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	 const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<Teuchos::ParameterList>& solverParams);

	//! Destructor
	virtual ~SalingerBordering();

	//! Set blocks in extended linear system
	/*!
	 * \param group [in] Underlying group representing J
	 * \param tpGroup [in] Turning point group representing the turning
	 *        point equations.  
	 * \param nullVector [in] Vector representing v
	 * \param JnVector [in] Vector representing Jv
	 * \param dfdp [in] Vector representing df/dp
	 * \param dJndp [in] Vector representing d(Jv)/dp
	 */
	virtual void setBlocks(
	  const Teuchos::RefCountPtr<LOCA::TurningPoint::MooreSpence::AbstractGroup>& group,
	  const Teuchos::RefCountPtr<LOCA::TurningPoint::MooreSpence::ExtendedGroup>& tpGroup,
	  const Teuchos::RefCountPtr<const NOX::Abstract::Vector>& nullVector,
	  const Teuchos::RefCountPtr<const NOX::Abstract::Vector>& JnVector,
	  const Teuchos::RefCountPtr<const NOX::Abstract::Vector>& dfdp,
	  const Teuchos::RefCountPtr<const NOX::Abstract::Vector>& dJndp);

	/*! 
	 * \brief Solves the extended system as defined above
	 */
	/*!
	 * The \em params argument is the linear solver parameters. 
	 */
	virtual NOX::Abstract::Group::ReturnType 
	solve(
	   Teuchos::ParameterList& params,
	   const LOCA::TurningPoint::MooreSpence::ExtendedMultiVector& input,
	   LOCA::TurningPoint::MooreSpence::ExtendedMultiVector& result) const;

      protected:

	//! Solves equations with contiguous arguments
	NOX::Abstract::Group::ReturnType 
	solveContiguous(
	          Teuchos::ParameterList& params,
		  const NOX::Abstract::MultiVector& input_x,
		  const NOX::Abstract::MultiVector& input_null,
	          const NOX::Abstract::MultiVector::DenseMatrix& input_param,
		  NOX::Abstract::MultiVector& result_x,
		  NOX::Abstract::MultiVector& result_null,
	          NOX::Abstract::MultiVector::DenseMatrix& result_param) const;

      private:

	//! Private to prohibit copying
	SalingerBordering(const SalingerBordering&);
	
	//! Private to prohibit copying
	SalingerBordering& operator = (const SalingerBordering&);

      protected:

	//! Global data object
	Teuchos::RefCountPtr<LOCA::GlobalData> globalData;
	
	//! Solver parameters
	Teuchos::RefCountPtr<Teuchos::ParameterList> solverParams;

	//! Underlying group
	Teuchos::RefCountPtr<LOCA::TurningPoint::MooreSpence::AbstractGroup> group;

	//! Turning point group
	Teuchos::RefCountPtr<LOCA::TurningPoint::MooreSpence::ExtendedGroup> tpGroup;

	//! Null vector
	Teuchos::RefCountPtr<const NOX::Abstract::Vector> nullVector;

	//! Jacobian times null vector
	Teuchos::RefCountPtr<const NOX::Abstract::Vector> JnVector;

	//! df/dp
	Teuchos::RefCountPtr<const NOX::Abstract::Vector> dfdp;

	//! d(Jn)/dp
	Teuchos::RefCountPtr<const NOX::Abstract::Vector> dJndp;

      }; // Class SalingerBordering
      
    } // Namespace MooreSpence

  } // Namespace TurningPoint
    
} // Namespace LOCA

#endif // LOCA_TURNINGPOINT_MOORESPENCE_SALINGERBORDERING_H
