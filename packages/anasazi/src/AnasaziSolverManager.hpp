// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
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
