// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_ITERATION_HPP
#define BELOS_ITERATION_HPP

/*! \file BelosIteration.hpp
    \brief Pure virtual base class which describes the basic interface to the linear solver iteration.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"

#include "Teuchos_Array.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"


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

  //! Destructor.
  virtual ~Iteration() {};
  //@}


  //! @name Solver methods
  //@{ 
  
  /*! \brief This method performs linear solver iterations until the status test
    indicates the need to stop or an error occurs (in which case, an std::exception is thrown).
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
  
  //! \brief Reset the iteration count to iter.
  virtual void resetNumIters( int iter = 0 ) = 0;

  //! Get the residuals native to the solver.
  //! \return A multivector with blockSize vectors containing the native residuals, else the native residual norm is returned.
  virtual Teuchos::RCP<const MV> getNativeResiduals( std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> *norms ) const = 0;

  //! Get the current update to the linear system.
  /*! \note Some solvers, like GMRES, do not compute updates to the solution every iteration.
            This method forces the computation of the current update.
  */
  virtual Teuchos::RCP<MV> getCurrentUpdate() const = 0;

  //@}


  
    //! @name Accessor methods
  //@{ 

  //! Get a constant reference to the linear problem.
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
