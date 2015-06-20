
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

  /// \brief Tolerance achieved by the last \c solve() invocation.
  ///
  /// This is the maximum over all right-hand sides' achieved
  /// convergence tolerances, and is set whether or not the solve
  /// actually managed to achieve the desired convergence tolerance.
  ///
  /// The default implementation throws std::runtime_error.  This is
  /// in case the idea of a single convergence tolerance doesn't make
  /// sense for some solvers.  It also serves as a gradual upgrade
  /// path (since this method is a later addition to the \c
  /// SolverManager interface).
  virtual typename Teuchos::ScalarTraits<ScalarType>::magnitudeType achievedTol() const {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "achievedTol() not implemented");
  }

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

  /// \brief Set the parameters to use when solving the linear problem.
  ///
  /// \param params [in/out] List of parameters to use when solving
  ///   the linear problem.  This list will be modified as necessary
  ///   to include default parameters that need not be provided.  If
  ///   params is null, then this method uses default parameters.
  ///
  /// \note The ParameterList returned by \c getValidParameters() has
  ///   all the parameters that the solver understands, possibly
  ///   including human-readable documentation and validators.
  virtual void setParameters( const Teuchos::RCP<Teuchos::ParameterList> &params ) = 0;

  //! Set user-defined convergence status test.
  virtual void setUserConvStatusTest(
    const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &userConvStatusTest
    )
    {
      TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Error, the function setUserConvStatusTest() has not been"
        << " overridden for the class" << this->description() << " yet!");
    }

  //! Set user-defined debug status test.
  virtual void setDebugStatusTest(
    const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &debugStatusTest
    )
    {
      TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Error, the function setDebugStatusTest() has not been"
        << " overridden for the class" << this->description() << " yet!");
    }

  //@}

  //! @name Reset methods
  //@{

  /// \brief Reset the solver manager.
  ///
  /// Reset the solver manager in a way specified by the \c
  /// ResetType parameter.  This informs the solver manager that the
  /// solver should prepare for the next call to solve by resetting
  /// certain elements of the iterative solver strategy.
  virtual void reset( const ResetType type ) = 0;
  //@}

  //! @name Solver application methods
  //@{

  /// \brief Iterate until the status test tells us to stop.
  //
  /// This method performs possibly repeated calls to the underlying
  /// linear solver's iterate() routine, until the problem has been
  /// solved (as decided by the solver manager via the status
  /// test(s)), or the solver manager decides to quit.
  ///
  /// \return A \c Belos::ReturnType enum specifying:
  ///   - Belos::Converged: the linear problem was solved to the
  ///     specification required by the solver manager.
  ///   - Belos::Unconverged: the linear problem was not solved to the
  ///     specification desired by the solver manager.
  virtual ReturnType solve() = 0;
  //@}

};


namespace Details {

  /// \class RealSolverManager
  /// \brief Base class for Belos::SolverManager subclasses which
  ///   normally can only compile for real ScalarType.
  ///
  /// It can be a lot of trouble to make a solver implementation
  /// compile correctly for both real and complex ScalarType.  That's
  /// why this class exists.  If a Belos solver inherits from this
  /// class, that indicates that the solver's implementation can only
  /// compile if ScalarType is real (not complex).  If you attempt to
  /// invoke the solver's constructor when <tt>ScalarType</tt> is
  /// complex, the constructor will throw an exception.
  ///
  /// The point of this class is to ensure that Belos::SolverFactory
  /// will always compile, even if some Belos solvers' implementations
  /// cannot compile when ScalarType is complex.  See GCRODRSolMgr for
  /// an example of how to use this class to avoid compilation of code
  /// that only works for real ScalarType.
  template<class ScalarType,
           class MV,
           class OP,
           const bool isComplex = Teuchos::ScalarTraits<ScalarType>::isComplex>
  class RealSolverManager;

  // Specialization for isComplex = true adds nothing to SolverManager.
  template<class ScalarType, class MV, class OP>
  class RealSolverManager<ScalarType, MV, OP, false> :
    public SolverManager<ScalarType, MV, OP> {
  public:
    RealSolverManager () {}
    virtual ~RealSolverManager () {}
  };

  // Specialization for isComplex = true (ScalarType is complex) adds
  // a constructor that always throws std::logic_error.  Subclasses
  // must always call the base class constructor.
  //
  // The complex version (isComplex = true) needs to implement all the
  // pure virtual methods in SolverManager, even though they can never
  // actually be called, since the constructor throws.
  template<class ScalarType, class MV, class OP>
  class RealSolverManager<ScalarType, MV, OP, true> :
    public SolverManager<ScalarType, MV, OP> {
  public:
    RealSolverManager () {
      TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
        "This solver is not implemented for complex ScalarType." );
    }
    virtual ~RealSolverManager () {}

    virtual const LinearProblem<ScalarType,MV,OP>& getProblem() const {
      TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
        "This solver is not implemented for complex ScalarType." );
    }
    virtual Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const {
      TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
        "This solver is not implemented for complex ScalarType." );
    }
    virtual Teuchos::RCP<const Teuchos::ParameterList> getCurrentParameters() const {
      TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
        "This solver is not implemented for complex ScalarType." );
    }
    virtual int getNumIters() const {
      TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
        "This solver is not implemented for complex ScalarType." );
    }
    virtual bool isLOADetected() const {
      TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
        "This solver is not implemented for complex ScalarType." );
    }
    virtual void setProblem (const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem) {
      TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
        "This solver is not implemented for complex ScalarType." );
    }
    virtual void setParameters (const Teuchos::RCP<Teuchos::ParameterList>& params) {
      TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
        "This solver is not implemented for complex ScalarType." );
    }
    virtual void reset (const ResetType type) {
      TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
        "This solver is not implemented for complex ScalarType." );
    }
    virtual ReturnType solve () {
      TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
        "This solver is not implemented for complex ScalarType." );
    }
  };

} // namespace Details


} // End Belos namespace

#endif /* BELOS_SOLVERMANAGER_HPP */
