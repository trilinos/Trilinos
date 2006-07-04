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

#ifndef LOCA_PITCHFORK_MOORESPENCE_SOLVERSTRATEGY_H
#define LOCA_PITCHFORK_MOORESPENCE_SOLVERSTRATEGY_H

#include "Teuchos_RefCountPtr.hpp" 

#include "NOX_Abstract_Group.H"  // for ReturnType

// forward declarations
namespace LOCA {
  namespace Pitchfork {
    namespace MooreSpence {
      class AbstractGroup;
      class ExtendedGroup;
      class ExtendedMultiVector;
    }
  }
}

namespace LOCA {

  namespace Pitchfork {

    namespace MooreSpence {
 
      /*! 
       * \brief Abstract strategy for solving the Moore-Spence pitchfork
       * equations.
       */
      /*!
       * This class provides an abstract interface for solver strategies to
       * solve the Moore-Spence pitchfork Newton system:
       * \f[
       *     \begin{bmatrix}
       *          J   &    0   & \psi &  f_p   \\
       *       (Jv)_x &    J   &   0  & (Jv)_p \\
       *       \psi^T &    0   &   0  &   0    \\
       *          0   & \phi^T &   0  &   0
       *     \end{bmatrix}
       *     \begin{bmatrix}
       *        X \\
       *        Y \\
       *        w \\
       *        z
       *     \end{bmatrix} = 
       *     \begin{bmatrix}
       *        F \\
       *        G \\
       *        s \\
       *        h
       *     \end{bmatrix}.
       *  \f]
       * After instantiating a solver Solvers (via
       * LOCA::Pitchfork::MooreSpence::SolverFactory), the linear
       * system is set up by setBlocks() and can then be solved by solve().
       */
      class SolverStrategy {

      public:

	//! Constructor
	SolverStrategy() {}

	//! Destructor
	virtual ~SolverStrategy() {}

	//! Set blocks in extended linear system
	/*!
	 * \param group [in] Underlying group representing J
	 * \param pfGroup [in] Pitchfork group representing the pitchfork 
	 *                     equations.  
	 * \param asymMultiVector [in] Multivector representing the 
	 *                             asymmetric vector
	 * \param nullVector [in] Vector representing v
	 * \param JnVector [in] Vector representing Jv
	 * \param dfdp [in] Vector representing df/dp
	 * \param dJndp [in] Vector representing d(Jv)/dp
	 */
	virtual void setBlocks(
	  const Teuchos::RefCountPtr<LOCA::Pitchfork::MooreSpence::AbstractGroup>& group,
	  const Teuchos::RefCountPtr<LOCA::Pitchfork::MooreSpence::ExtendedGroup>& pfGroup,
	  const Teuchos::RefCountPtr<const NOX::Abstract::MultiVector>& asymMultiVector,
	  const Teuchos::RefCountPtr<const NOX::Abstract::Vector>& nullVector,
	  const Teuchos::RefCountPtr<const NOX::Abstract::Vector>& JnVector,
	  const Teuchos::RefCountPtr<const NOX::Abstract::Vector>& dfdp,
	  const Teuchos::RefCountPtr<const NOX::Abstract::Vector>& dJndp) = 0;

	/*! 
	 * \brief Solves the extended system as defined above
	 */
	/*!
	 * The \em params argument is the linear solver parameters. 
	 */
	virtual NOX::Abstract::Group::ReturnType 
	solve(
	  Teuchos::ParameterList& params,
	  const LOCA::Pitchfork::MooreSpence::ExtendedMultiVector& input,
          LOCA::Pitchfork::MooreSpence::ExtendedMultiVector& result) const = 0;

      private:

	//! Private to prohibit copying
	SolverStrategy(const SolverStrategy&);
	
	//! Private to prohibit copying
	SolverStrategy& operator = (const SolverStrategy&);

      }; // Class SolverStrategy
      
    } // Namespace MooreSpence

  } // Namespace TurningPont
    
} // Namespace LOCA

#endif // LOCA_PITCHFORK_MOORESPENCE_SOLVERSTRATEGY_H
