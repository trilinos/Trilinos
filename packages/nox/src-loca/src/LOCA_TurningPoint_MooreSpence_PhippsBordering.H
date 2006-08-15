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

#ifndef LOCA_TURNINGPOINT_MOORESPENCE_PHIPPSBORDERING_H
#define LOCA_TURNINGPOINT_MOORESPENCE_PHIPPSBORDERING_H

#include "LOCA_TurningPoint_MooreSpence_SolverStrategy.H"
#include "NOX_Abstract_MultiVector.H"

// forward declarations
namespace LOCA {
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }
  namespace BorderedSolver {
    class AbstractStrategy;
  }
}

namespace LOCA {

  namespace TurningPoint {

    namespace MooreSpence {
 
      /*!
       * \brief Moore-Spence turning point solver strategy based on
       * "Phipps" bordering which is the 5-solve modified turning
       * point bordering algorithm that uses bordered linear solves.
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
       *  via the following modified block elimination scheme:
       *  \f[
       *    \begin{split}
       *      \begin{bmatrix} 
       *         J   & u \\
       *         v^T & 0
       *      \end{bmatrix}
       *      \begin{bmatrix}
       *         A & B \\
       *         a & b
       *      \end{bmatrix} &= 
       *      \begin{bmatrix}
       *         F & f_p \\
       *         0 & 0
       *      \end{bmatrix} \\
       *      \begin{bmatrix} 
       *         J   & u \\
       *         v^T & 0
       *      \end{bmatrix}
       *      \begin{bmatrix}
       *         C & D & E\\
       *         c & d & e
       *      \end{bmatrix} &= 
       *      \begin{bmatrix}
       *         (Jv)_x A - G & (Jv)_x B - (Jv)_p & (Jv)_x v \\
       *                0    &          0         &     0    
       *      \end{bmatrix} \\
       *      \begin{bmatrix}
       *         \sigma    &    0     &  b      \\
       *            e      & \sigma   & -d      \\
       *         -\phi^T E & \phi^T v & \phi^T D
       *      \end{bmatrix}
       *      \begin{bmatrix}
       *         \alpha \\
       *         \beta \\
       *         z
       *      \end{bmatrix} &=
       *      \begin{bmatrix}
       *         a \\
       *         c \\
       *         h + \phi^T C
       *      \end{bmatrix} \\
       *       X &= A - B z + v \alpha \\
       *       Y &= -C + d z - E \alpha + v \beta
       *    \end{split}
       *   \f]
       * where \f$s = \|J v\|\f$ and \f$u = J v/s\f$.  Each bordered solve
       * is implemented by a LOCA::BorderedSolver::AbstractStrategy strategy
       * object.
       */
      class PhippsBordering : 
	public LOCA::TurningPoint::MooreSpence::SolverStrategy {

      public:

	//! Constructor. 
	/*!
	 * \param global_data [in] Global data object
	 * \param topParams [in] Parsed top-level parameter list 
	 * \param solverParams [in] Bordered solver parameters.  Instantiates 
	 *        a bordered solver for solving the bordeded systems 
	 *        described above.  See LOCA::BorderedSolver::Factory for a
	 *        description of available solvers.
	 */
	PhippsBordering(
	 const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	 const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<Teuchos::ParameterList>& solverParams);

	//! Destructor
	virtual ~PhippsBordering();

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
	PhippsBordering(const PhippsBordering&);
	
	//! Private to prohibit copying
	PhippsBordering& operator = (const PhippsBordering&);

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

	//! Bordered solver for solving (n+1)x(n+1) sets of equations
	Teuchos::RefCountPtr<LOCA::BorderedSolver::AbstractStrategy> borderedSolver;

	//! Null vector as a multivector
	Teuchos::RefCountPtr<NOX::Abstract::MultiVector> nullMultiVector;

	//! Jacobian times null vector as a multivector
	Teuchos::RefCountPtr<NOX::Abstract::MultiVector> JnMultiVector;

	//! Norm of JnVector
	double s;

      }; // Class PhippsBordering
      
    } // Namespace MooreSpence

  } // Namespace TurningPoint
    
} // Namespace LOCA

#endif // LOCA_TURNINGPOINT_MOORESPENCE_PHIPPSBORDERING_H
