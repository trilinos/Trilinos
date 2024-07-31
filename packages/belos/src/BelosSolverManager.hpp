// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_SOLVERMANAGER_HPP
#define BELOS_SOLVERMANAGER_HPP

/*! \file BelosSolverManager.hpp
    \brief Pure virtual base class which describes the basic interface for a solver manager.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosStatusTestCombo.hpp"

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

  /// \brief clone the solver manager.
  ///
  /// Implements the DII inversion and injection pattern
  virtual Teuchos::RCP<SolverManager<ScalarType, MV, OP> > clone () const = 0;
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
    const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &/* userConvStatusTest */,
    const typename StatusTestCombo<ScalarType,MV,OP>::ComboType &/* comboType */ =
        StatusTestCombo<ScalarType,MV,OP>::SEQ
    )
    {
      TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Error, the function setUserConvStatusTest() has not been"
        << " overridden for the class" << this->description() << " yet!");
    }

  //! Set user-defined debug status test.
  virtual void setDebugStatusTest(
    const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &/* debugStatusTest */
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
      // Do not throw on constructor. The DII system registers all class types
      // and must construct even if the class will not be usable.
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


  /// \class LapackSupportsScalar
  /// \brief Type traits class that says whether Teuchos::LAPACK
  ///   has a valid implementation for the given ScalarType.
  template<class ScalarType>
  class LapackSupportsScalar {
  public:
    const static bool value = false;
  };

  template<>
  class LapackSupportsScalar<float> {
  public:
    const static bool value = true;
  };

  template<>
  class LapackSupportsScalar<double> {
  public:
    const static bool value = true;
  };
  
#ifdef HAVE_TEUCHOS_LONG_DOUBLE
  template<>
  class LapackSupportsScalar<long double> {
  public:
    const static bool value = true;
  };
#endif


#ifdef HAVE_TEUCHOS_COMPLEX
  template<>
  class LapackSupportsScalar<std::complex<float> > {
  public:
    const static bool value = true;
  };

  template<>
  class LapackSupportsScalar<std::complex<double> > {
  public:
    const static bool value = true;
  };
#endif // HAVE_TEUCHOS_COMPLEX

  /// \class SolverManagerRequiresLapack
  /// \brief Base class for Belos::SolverManager subclasses which
  ///   normally can only compile with ScalarType types for which
  ///   Teuchos::LAPACK has a valid implementation.
  template<class ScalarType,
           class MV,
           class OP,
           const bool lapackSupportsScalarType =
           Belos::Details::LapackSupportsScalar<ScalarType>::value>
  class SolverManagerRequiresLapack;

  /// \brief Specialization for ScalarType types for which
  ///   Teuchos::LAPACK has a valid implementation.
  ///
  /// This specialization adds nothing to SolverManager.
  template<class ScalarType, class MV, class OP>
  class SolverManagerRequiresLapack<ScalarType, MV, OP, true> :
    public SolverManager<ScalarType, MV, OP> {
  public:
    SolverManagerRequiresLapack () {}
    virtual ~SolverManagerRequiresLapack () {}
  };

  /// \brief Specialization for ScalarType types for which
  ///   Teuchos::LAPACK does NOT have a valid implementation.
  ///
  /// This is a stub specialization whose constructor always throws
  /// std::logic_error.  Subclasses must always call the base class
  /// constructor.
  template<class ScalarType, class MV, class OP>
  class SolverManagerRequiresLapack<ScalarType, MV, OP, false> :
    public SolverManager<ScalarType, MV, OP> {
  public:
    SolverManagerRequiresLapack () {
      TEUCHOS_TEST_FOR_EXCEPTION
        (true, std::logic_error, "This solver is not implemented for ScalarType"
         " types for which Teuchos::LAPACK does not have a valid implementation.  "
         "ScalarType = " << Teuchos::TypeNameTraits<ScalarType>::name () << ".");
    }
    virtual ~SolverManagerRequiresLapack () {}

    virtual const LinearProblem<ScalarType,MV,OP>& getProblem() const {
      TEUCHOS_TEST_FOR_EXCEPTION
        (true, std::logic_error, "This solver is not implemented for ScalarType"
         " types for which Teuchos::LAPACK does not have a valid implementation.  "
         "ScalarType = " << Teuchos::TypeNameTraits<ScalarType>::name () << ".");
    }
    virtual Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const {
      TEUCHOS_TEST_FOR_EXCEPTION
        (true, std::logic_error, "This solver is not implemented for ScalarType"
         " types for which Teuchos::LAPACK does not have a valid implementation.  "
         "ScalarType = " << Teuchos::TypeNameTraits<ScalarType>::name () << ".");
    }
    virtual Teuchos::RCP<const Teuchos::ParameterList> getCurrentParameters() const {
      TEUCHOS_TEST_FOR_EXCEPTION
        (true, std::logic_error, "This solver is not implemented for ScalarType"
         " types for which Teuchos::LAPACK does not have a valid implementation.  "
         "ScalarType = " << Teuchos::TypeNameTraits<ScalarType>::name () << ".");
    }
    virtual int getNumIters() const {
      TEUCHOS_TEST_FOR_EXCEPTION
        (true, std::logic_error, "This solver is not implemented for ScalarType"
         " types for which Teuchos::LAPACK does not have a valid implementation.  "
         "ScalarType = " << Teuchos::TypeNameTraits<ScalarType>::name () << ".");
    }
    virtual bool isLOADetected() const {
      TEUCHOS_TEST_FOR_EXCEPTION
        (true, std::logic_error, "This solver is not implemented for ScalarType"
         " types for which Teuchos::LAPACK does not have a valid implementation.  "
         "ScalarType = " << Teuchos::TypeNameTraits<ScalarType>::name () << ".");
    }
    virtual void setProblem (const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem) {
      TEUCHOS_TEST_FOR_EXCEPTION
        (true, std::logic_error, "This solver is not implemented for ScalarType"
         " types for which Teuchos::LAPACK does not have a valid implementation.  "
         "ScalarType = " << Teuchos::TypeNameTraits<ScalarType>::name () << ".");
    }
    virtual void setParameters (const Teuchos::RCP<Teuchos::ParameterList>& params) {
      TEUCHOS_TEST_FOR_EXCEPTION
        (true, std::logic_error, "This solver is not implemented for ScalarType"
         " types for which Teuchos::LAPACK does not have a valid implementation.  "
         "ScalarType = " << Teuchos::TypeNameTraits<ScalarType>::name () << ".");
    }
    virtual void reset (const ResetType type) {
      TEUCHOS_TEST_FOR_EXCEPTION
        (true, std::logic_error, "This solver is not implemented for ScalarType"
         " types for which Teuchos::LAPACK does not have a valid implementation.  "
         "ScalarType = " << Teuchos::TypeNameTraits<ScalarType>::name () << ".");
    }
    virtual ReturnType solve () {
      TEUCHOS_TEST_FOR_EXCEPTION
        (true, std::logic_error, "This solver is not implemented for ScalarType"
         " types for which Teuchos::LAPACK does not have a valid implementation.  "
         "ScalarType = " << Teuchos::TypeNameTraits<ScalarType>::name () << ".");
    }
  };

  /// \class SolverManagerRequiresRealLapack
  /// \brief Base class for Belos::SolverManager subclasses which
  ///   normally can only compile with <i>real</i> ScalarType types
  ///   for which Teuchos::LAPACK has a valid implementation.
  template<class ScalarType,
           class MV,
           class OP,
           const bool supportsScalarType =
             Belos::Details::LapackSupportsScalar<ScalarType>::value &&
             ! Teuchos::ScalarTraits<ScalarType>::isComplex>
  class SolverManagerRequiresRealLapack;

  /// \brief Non-stub specialization for real ScalarType types for
  ///   which Teuchos::LAPACK has a valid implementation.
  ///
  /// This specialization adds nothing to SolverManager.  Thus, the
  /// SolverManager subclass that has the actual specific solver
  /// implementation gets to implement any virtual methods of
  /// SolverManager.
  template<class ScalarType, class MV, class OP>
  class SolverManagerRequiresRealLapack<ScalarType, MV, OP, true> :
    public SolverManager<ScalarType, MV, OP> {
  public:
    SolverManagerRequiresRealLapack () {}
    virtual ~SolverManagerRequiresRealLapack () {}
  };

  /// \brief Stub specialization for ScalarType types which are NOT
  ///   real, or for which Teuchos::LAPACK does NOT have a valid
  ///   implementation.
  ///
  /// This is a stub specialization whose constructor always throws
  /// std::logic_error.  Subclasses must always call the base class
  /// constructor.
  template<class ScalarType, class MV, class OP>
  class SolverManagerRequiresRealLapack<ScalarType, MV, OP, false> :
    public SolverManager<ScalarType, MV, OP> {
  public:
    SolverManagerRequiresRealLapack () {
      // Do not throw on constructor. The DII system registers all class types
      // and must construct even if the class will not be usable.
    }
    virtual ~SolverManagerRequiresRealLapack () {}

    virtual const LinearProblem<ScalarType,MV,OP>& getProblem() const {
      TEUCHOS_TEST_FOR_EXCEPTION
        (true, std::logic_error, "This solver is not implemented for complex "
         "ScalarType types, or for ScalarType types for which Teuchos::LAPACK "
         "does not have a valid implementation."
         "ScalarType = " << Teuchos::TypeNameTraits<ScalarType>::name () << ".");
    }
    virtual Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const {
      TEUCHOS_TEST_FOR_EXCEPTION
        (true, std::logic_error, "This solver is not implemented for complex "
         "ScalarType types, or for ScalarType types for which Teuchos::LAPACK "
         "does not have a valid implementation."
         "ScalarType = " << Teuchos::TypeNameTraits<ScalarType>::name () << ".");
    }
    virtual Teuchos::RCP<const Teuchos::ParameterList> getCurrentParameters() const {
      TEUCHOS_TEST_FOR_EXCEPTION
        (true, std::logic_error, "This solver is not implemented for complex "
         "ScalarType types, or for ScalarType types for which Teuchos::LAPACK "
         "does not have a valid implementation."
         "ScalarType = " << Teuchos::TypeNameTraits<ScalarType>::name () << ".");
    }
    virtual int getNumIters() const {
      TEUCHOS_TEST_FOR_EXCEPTION
        (true, std::logic_error, "This solver is not implemented for complex "
         "ScalarType types, or for ScalarType types for which Teuchos::LAPACK "
         "does not have a valid implementation."
         "ScalarType = " << Teuchos::TypeNameTraits<ScalarType>::name () << ".");
    }
    virtual bool isLOADetected() const {
      TEUCHOS_TEST_FOR_EXCEPTION
        (true, std::logic_error, "This solver is not implemented for complex "
         "ScalarType types, or for ScalarType types for which Teuchos::LAPACK "
         "does not have a valid implementation."
         "ScalarType = " << Teuchos::TypeNameTraits<ScalarType>::name () << ".");
    }
    virtual void
    setProblem (const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> >& /* problem */) {
      TEUCHOS_TEST_FOR_EXCEPTION
        (true, std::logic_error, "This solver is not implemented for complex "
         "ScalarType types, or for ScalarType types for which Teuchos::LAPACK "
         "does not have a valid implementation."
         "ScalarType = " << Teuchos::TypeNameTraits<ScalarType>::name () << ".");
    }
    virtual void setParameters (const Teuchos::RCP<Teuchos::ParameterList>& params) {
      TEUCHOS_TEST_FOR_EXCEPTION
        (true, std::logic_error, "This solver is not implemented for complex "
         "ScalarType types, or for ScalarType types for which Teuchos::LAPACK "
         "does not have a valid implementation."
         "ScalarType = " << Teuchos::TypeNameTraits<ScalarType>::name () << ".");
    }
    virtual void reset (const ResetType type) {
      TEUCHOS_TEST_FOR_EXCEPTION
        (true, std::logic_error, "This solver is not implemented for complex "
         "ScalarType types, or for ScalarType types for which Teuchos::LAPACK "
         "does not have a valid implementation."
         "ScalarType = " << Teuchos::TypeNameTraits<ScalarType>::name () << ".");
    }
    virtual ReturnType solve () {
      TEUCHOS_TEST_FOR_EXCEPTION
        (true, std::logic_error, "This solver is not implemented for complex "
         "ScalarType types, or for ScalarType types for which Teuchos::LAPACK "
         "does not have a valid implementation."
         "ScalarType = " << Teuchos::TypeNameTraits<ScalarType>::name () << ".");
    }
  };

} // namespace Details
} // namespace Belos

#endif /* BELOS_SOLVERMANAGER_HPP */
