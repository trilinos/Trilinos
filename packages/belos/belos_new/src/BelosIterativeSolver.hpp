
// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef BELOS_ITERATIVE_SOLVER_HPP
#define BELOS_ITERATIVE_SOLVER_HPP

/*! \file BelosIterativeSolver.hpp
    \brief Pure virtual base class which describes the basic interface to the iterative solver.
*/

#include "BelosLinearProblemManager.hpp"
#include "BelosReturnType.hpp"
#include "BelosOperator.hpp"
#include "BelosMultiVec.hpp"

/*! \class Belos::IterativeSolver
  \brief The Belos::IterativeSolver is a templated virtual base class that defines the
	basic interface that any linear solver will support.

	The Belos::IterativeSolver class is responsible for providing the 
	current solver information to the Belos::StatusTest object.
*/

namespace Belos {

template<class TYPE>
class IterativeSolver {
    
  public:

  //@{ \name Constructor/Destructor.

  //! Default Constructor.
  IterativeSolver(void) {};

  //! Destructor.
  virtual ~IterativeSolver(void) {};
  //@}
  
  //@{ \name Accessor methods

  //! Get the current iteration count for this block of linear systems.
  virtual int GetNumIters() const = 0;

  //! Get the current restart count of the iteration method.
  /*! Some linear solvers can perform restarts (i.e. GMRES) to reduce memory
	and orthogonalization costs.  For other linear solvers that don't
	perform restarts (i.e. CG), this is not a valid stopping criteria.
  */
  virtual int GetNumRestarts() const = 0;

  //! Get the solvers native residuals for the current block of linear systems. 
  /*! This is not be the same as the true residual for most solvers. Somtimes the native
    residuals are not in multivector form, so the norm type is solver dependent.  
    If the true residual vectors are required, then call GetTrueResidVecs().
    \note The norm requested for this operation may not be supported by the solver.
    If it is not, the return type will be undefined.
  */
  virtual ReturnType GetNativeResidNorms(TYPE* normvec, NormType norm_type) const = 0;

  //! Get the actual residual vectors for the current block of linear systems.
  /*! This may force the solver to compute a current residual for its linear
  	systems.  Some linear solvers don't constantly update their 
	residuals and solutions (i.e. GMRES), so this may be an expensive request.
	<b>Using true residuals to determine convergence should be secondary
	to using the native residuals of the iterative linear solver.</b>
  */
  virtual MultiVec<TYPE>* GetCurrentSoln() = 0;

  /*! \brief Get a constant reference to the current linear problem, 
    	which may include a current solution.
  */
  virtual LinearProblemManager<TYPE>& GetLinearProblem() const = 0;

  //@}

};

} // end Belos namespace

#endif /* BELOS_ITERATIVE_SOLVER_HPP */
