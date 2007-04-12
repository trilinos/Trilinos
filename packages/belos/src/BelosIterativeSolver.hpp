
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

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_Describable.hpp"

using Teuchos::RefCountPtr;
using Teuchos::ParameterList;


/*! \class Belos::IterativeSolver
  \brief The Belos::IterativeSolver is a templated virtual base class that defines the
	basic interface that any linear solver will support.

	The Belos::IterativeSolver class is responsible for providing the 
	current solver information to the Belos::StatusTest object.
*/

namespace Belos {

template <class ScalarType>
class OutputManager;

template <class ScalarType, class MV, class OP>
class StatusTest;

template <class ScalarType, class MV, class OP>
class LinearProblem;

template <class ScalarType, class MV, class OP>
class IterativeSolver : virtual public Teuchos::Describable {
 
  typedef Teuchos::ScalarTraits<ScalarType> SCT;
  typedef typename SCT::magnitudeType MagnitudeType;
  
  public:

  /** \name Constructor/Destructor. */
  //@{

  //! Default Constructor.
  IterativeSolver(void) {};

  //! Destructor.
  virtual ~IterativeSolver(void) {};

  //@}
  
  /** \name Accessor methods */
  //@{

  //! Get the current iteration count for this block of linear systems.
  virtual int GetNumIters() const = 0;

  //! Get the current restart count of the iteration method.
  /*! Some linear solvers can perform restarts (i.e. GMRES) to reduce memory
	and orthogonalization costs.  For other linear solvers that don't
	perform restarts (i.e. CG), this is not a valid stopping criteria.
  */
  virtual int GetNumRestarts() const = 0;

  //! Get the solvers native residuals for the current block of linear systems. 
  /*! This is not be the same as the true residual for most solvers. Sometimes the native
    residuals are not in multivector form, so the norm type is solver dependent.  
    If the true residual is required, then call GetCurrentSoln().

    \note
    <ol>
      <li> If the native residual is in multivector form then a non-null pointer will be
      returned, else the normvec will be populated with the current residual norms. 
      <li> If the native residual is returned in multivector form, the memory is managed
      by the calling routine.
    </ol>
  */
  virtual RefCountPtr<const MV> GetNativeResiduals( std::vector<MagnitudeType> *normvec ) const = 0;

  //! Get the actual residual vectors for the current block of linear systems.
  /*! This may force the solver to compute a current residual for its linear
    systems.  Some linear solvers don't constantly update their 
    residuals and solutions (i.e. GMRES), so this may be an expensive request.
    <b>Using true residuals to determine convergence should be secondary
    to using the native residuals of the iterative linear solver.</b>

    \note The memory of the returned multivector is managed by the calling routine.
  */
  virtual RefCountPtr<MV> GetCurrentSoln() = 0;

  /*! \brief Get a constant reference to the current linear problem, 
    	which may include a current solution.
  */
  virtual RefCountPtr<LinearProblem<ScalarType,MV,OP> > GetLinearProblem() const = 0;

  virtual RefCountPtr<StatusTest<ScalarType,MV,OP> > GetStatusTest() const = 0;

  //@}

  /** \name Reset methods */
  //@{
  
  /*! \brief Reset the solver to its initialized state.
     This is not a required method for all solvers, but a method that should be
     implemented by any solver whos object maintains state information.  The reset
     method is used to reset a solver object so it can then be reused, i.e. in
     inner-outer iterations.  An optional parameter list can be passed in to change
     any solver parameters necessary.
     \note This method can ONLY be expected to reset the solver object's state, 
     the LinearProblem and StatusTest need to be reset separately.
  */
  virtual int Reset( const RefCountPtr<ParameterList>& pl = Teuchos::null,
                     const RefCountPtr<LinearProblem<ScalarType,MV,OP> > &lp = Teuchos::null,
                     const RefCountPtr<StatusTest<ScalarType,MV,OP> > &stest = Teuchos::null,
                     const RefCountPtr<OutputManager<ScalarType> > &om = Teuchos::null )
  { return 0; }
    
  //@}

  /** \name Solve method */
  //@{
 
  /*! \brief Use the iterative method prescribed by the solver to compute the solution
      to the linear problem.
   */
  virtual void Solve() = 0;

  //@}
};

} // end Belos namespace

#endif /* BELOS_ITERATIVE_SOLVER_HPP */
