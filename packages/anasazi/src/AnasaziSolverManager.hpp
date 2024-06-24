// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ANASAZI_SOLVERMANAGER_HPP
#define ANASAZI_SOLVERMANAGER_HPP

/*! \file AnasaziSolverManager.hpp
    \brief Pure virtual base class which describes the basic interface for a solver manager.
*/

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"
#include "AnasaziEigenproblem.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Time.hpp"

/*! \class Anasazi::SolverManager
  \brief The Anasazi::SolverManager is a templated virtual base class that defines the
	basic interface that any solver manager will support.
*/

namespace Anasazi {

template<class ScalarType, class MV, class OP>
class SolverManager {
    
  public:

    //!@name Constructors/Destructor 
  //@{ 

  //! Empty constructor.
  SolverManager() {};

  //! Destructor.
  virtual ~SolverManager() {};
  //@}
  
  //! @name Accessor methods
  //@{ 

  //! Return the eigenvalue problem.
  virtual const Eigenproblem<ScalarType,MV,OP>& getProblem() const = 0;

  //! Get the iteration count for the most recent call to \c solve().
  virtual int getNumIters() const = 0;

  //! Return the timers for this object. 
  virtual Teuchos::Array<Teuchos::RCP<Teuchos::Time> > getTimers() const 
  {
    Teuchos::Array<Teuchos::RCP<Teuchos::Time> > null_time;
    return null_time;
  }

  //@}

  //! @name Solver application methods
  //@{ 
    
  /*! \brief This method performs possibly repeated calls to the underlying eigensolver's iterate() routine
   * until the problem has been solved (as decided by the solver manager) or the solver manager decides to 
   * quit.
   *
   * \returns ::ReturnType specifying:
   *    - ::Converged: the eigenproblem was solved to the specification required by the solver manager.
   *    - ::Unconverged: the eigenproblem was not solved to the specification desired by the solver manager
  */
  virtual ReturnType solve() = 0;
  //@}
  
};

} // end Anasazi namespace

#endif /* ANASAZI_SOLVERMANAGER_HPP */
