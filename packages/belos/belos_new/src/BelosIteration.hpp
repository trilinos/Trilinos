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

#ifndef BELOS_ITERATION_HPP
#define BELOS_ITERATION_HPP

/*! \file BelosIteration.hpp
    \brief Pure virtual base class which describes the basic interface to the linear solver iteration.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Array.hpp"

namespace Belos {

template <class ScalarType, class MV, class OP>
class LinearProblem;

template <class ScalarType>
class OutputManager;

template <class ScalarType, class MV, class OP>
class StatusTest;

template <class ScalarType, class MV, class OP>
class MatOrthoManager;

template<class ScalarType, class MV, class OP>
class Iteration {

  public:

  //! @name Constructors/Destructor
  //@{ 

  //! Default Constructor.
  Iteration() {};

  //! Basic Constructor.
  /*! This constructor, implemented by all Belos iterations, takes an Belos::LinearProblem,
    Belos::OrthoManager, Belos::OutputManager, and Teuchos::ParameterList as input.  
    These four arguments are sufficient enough for constructing any Belos::Iteration object.
  */
  Iteration( const Teuchos::RefCountPtr<LinearProblem<ScalarType,MV,OP> > &problem, 
	     const Teuchos::RefCountPtr<OutputManager<ScalarType> > &printer,
	     const Teuchos::RefCountPtr<StatusTest<ScalarType,MV,OP> > &tester,
	     const Teuchos::RefCountPtr<MatOrthoManager<ScalarType,MV,OP> > &ortho,
	     Teuchos::ParameterList &params );

  //! Destructor.
  virtual ~Iteration() {};
  //@}


  //! @name Solver methods
  //@{ 
  
  /*! \brief This method performs linear solver iterations until the status test
    indicates the need to stop or an error occurs (in which case, an exception is thrown).
  */
  virtual void iterate() = 0;

  /*! \brief Initialize the solver with the initial vectors from the linear problem
   *  or random data.
   */
  virtual void initialize() = 0;

  //@}

  
  //! @name Status methods
  //@{ 

  //! \brief Get the current iteration count.
  virtual int getNumIters() const = 0;
  
  //! \brief Reset the iteration count.
  virtual void resetNumIters() = 0;

  //! Get the residuals native to the solver.
  //! \return A multivector with blockSize vectors containing the native residuals, else the native residual norm is returned.
  virtual Teuchos::RefCountPtr<const MV> getNativeResiduals( std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> *norms ) const = 0;

  //! Get the current update to the linear system.
  /*! \note Some solvers, like GMRES, do not compute updates to the solution every iteration.
            This method forces the computation of the current update.
  */
  virtual Teuchos::RefCountPtr<MV> getCurrentUpdate() const = 0;

  //! Get the dimension of the search subspace used to generate the current eigenvectors and eigenvalues.
  virtual int getCurSubspaceDim() const = 0;

  //! Get the maximum dimension allocated for the search subspace.
  virtual int getMaxSubspaceDim() const = 0;

  //@}


  
    //! @name Accessor methods
  //@{ 

  //! Get a constant reference to the eigenvalue problem.
  virtual const LinearProblem<ScalarType,MV,OP>& getProblem() const = 0;

  //! Get the blocksize to be used by the iterative solver in solving this linear problem.
  virtual int getBlockSize() const = 0;
  
  //! \brief Set the blocksize to be used by the iterative solver in solving this linear problem.
  virtual void setBlockSize(int blockSize) = 0;

  //! States whether the solver has been initialized or not.
  virtual bool isInitialized() = 0;

  //@}

};

} // end Belos namespace

#endif /* BELOS_ITERATION_HPP */
