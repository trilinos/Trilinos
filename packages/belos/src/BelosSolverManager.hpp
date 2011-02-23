
//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// ************************************************************************
//@HEADER

#ifndef BELOS_SOLVERMANAGER_HPP
#define BELOS_SOLVERMANAGER_HPP

/*! \file BelosSolverManager.hpp
    \brief Pure virtual base class which describes the basic interface for a solver manager.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"
#include "BelosLinearProblem.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Describable.hpp"

/*! \class Belos::SolverManager
  \brief The Belos::SolverManager is a templated virtual base class that defines the
	basic interface that any solver manager will support.
*/

namespace Belos {


template <class ScalarType, class MV, class OP>
class StatusTest;


template<class ScalarType, class MV, class OP>
class SolverManager : virtual public Teuchos::Describable {
    
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

  //! Return a reference to the linear problem being solved by this solver manager.
  virtual const LinearProblem<ScalarType,MV,OP>& getProblem() const = 0;

  //! Return the valid parameters for this solver manager.
  virtual Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const = 0;

  //! Return the current parameters being used for this solver manager.
  virtual Teuchos::RCP<const Teuchos::ParameterList> getCurrentParameters() const = 0;

  //! Get the iteration count for the most recent call to \c solve().
  virtual int getNumIters() const = 0;

  /*! \brief Returns whether a loss of accuracy was detected in the solver. 
   *  \note This method is normally applicable to GMRES-type solvers.
  */
  virtual bool isLOADetected() const = 0;
 
  //@}

  //! @name Set methods
  //@{

  //! Set the linear problem that needs to be solved. 
  virtual void setProblem( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem ) = 0;

  //! Set the parameters the solver manager should use to solve the linear problem.
  virtual void setParameters( const Teuchos::RCP<Teuchos::ParameterList> &params ) = 0;

  //! Set user-defined convergence status test.
  virtual void setUserConvStatusTest(
    const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &userConvStatusTest
    )
    {
      TEST_FOR_EXCEPT_MSG(true, "Error, the function setUserConvStatusTest() has not been"
        << " overridden for the class" << this->description() << " yet!");
    }

  //@}

  //! @name Reset methods
  //@{
  /*! \brief Performs a reset of the solver manager specified by the \c ResetType.  This informs the
   *  solver manager that the solver should prepare for the next call to solve by resetting certain elements
   *  of the iterative solver strategy.
  */ 
  virtual void reset( const ResetType type ) = 0;
  //@}

  //! @name Solver application methods
  //@{ 
    
  /*! \brief This method performs possibly repeated calls to the underlying linear solver's iterate() routine
   * until the problem has been solved (as decided by the solver manager) or the solver manager decides to 
   * quit.
   *
   * \returns ::ReturnType specifying:
   *    - ::Converged: the linear problem was solved to the specification required by the solver manager.
   *    - ::Unconverged: the linear problem was not solved to the specification desired by the solver manager
  */
  virtual ReturnType solve() = 0;
  //@}
  
};

} // end Belos namespace

#endif /* BELOS_SOLVERMANAGER_HPP */
