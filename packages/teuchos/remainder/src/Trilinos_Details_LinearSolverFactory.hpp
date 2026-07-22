// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TRILINOS_DETAILS_LINEARSOLVER_FACTORY_HPP
#define TRILINOS_DETAILS_LINEARSOLVER_FACTORY_HPP

/// \file Trilinos_Details_LinearSolverFactory.hpp
/// \brief Declaration and definition of linear solver factory, and
///   "factory of factories"
///
/// \warning This header file is NOT currently part of the public
///   interface of Trilinos.  It or its contents may change or
///   disappear at any time.
///
/// Tpetra::Details::getLinearSolver,
/// Tpetra::Details::registerLinearSolverFactory, and
/// Tpetra::Details::LinearSolverFactory implement the Dependency
/// Inversion and Injection (DII) pattern, as applied to "linear
/// solvers."  A linear solver solves or helps solve linear system(s)
/// AX=B.  Examples include sparse direct solvers, iterative solvers,
/// and preconditioners for iterative solvers.
///
/// DII naturally admits hierarchical run-time options, as in e.g.,
/// Teuchos::ParameterList.  This lets solvers create inner solvers in
/// an arbitrarily nested way, following the arbitrary nesting of the
/// Teuchos::ParameterList.
///
/// DII works well when a ParameterList can express all the data that
/// a solver might need.  However, some solvers need or may benefit
/// from additional data.  For example, algebraic multigrid can use
/// mesh coordinates, and a sparse factorization can use an initial
/// permutation.  Such data do not fit naturally in a
/// Teuchos::ParameterList.
///
/// \note To developers: The LinearSolver and LinearSolverFactory
///   interfaces, and the LinearSolverFactoryRepository interface and
///   implementation must live in the bottom-most (most upstream)
///   package from all solvers that depend on it.  Solver defines an
///   interface for a solver, and LinearSolverFactory defines an
///   interface for a "factory" that knows how to create solvers.
///   Each solver package defines its own solvers, and its own factory
///   that knows how to create all the solvers in a given package.

#include "Teuchos_RCP.hpp" // includes Teuchos_ConfigDefs.hpp
#include "TeuchosRemainder_config.h"
#include <map>
#ifdef HAVE_TEUCHOSCORE_CXX11
#  include <memory> // std::shared_ptr
#endif // HAVE_TEUCHOSCORE_CXX11
#include <stdexcept>
#include <sstream>
#include <string>


// Attempted fix for Bug 6392: declare all packages'
// LinearSolverFactory registration functions here, with weak linkage.
// This works whether or not the packages in question are actually
// enabled.  In createPackageNames() below, actually call these
// functions if they are linked in.  We only need to do this if
// building with static libraries; if building with dynamic libraries,
// each package takes care of this on its own.
//
// I wrote "attempted" because it DOESN'T WORK.  It doesn't matter
// whether these or their uses are in the .cpp or .hpp file, or
// whether they are in a regular function that gets compiled or a
// templated function that might not.
#if ! defined(HAVE_TEUCHOS_DYNAMIC_LIBS) && defined(HAVE_TEUCHOS_CXX_ATTRIBUTE_WEAK)
// FIXME (mfh 21 Aug 2015) NONE of the commented-out things work.

// namespace Amesos2 {
// namespace Details {
//   extern void __attribute__((weak)) registerLinearSolverFactory ();
// } // namespace Details
// } // namespace Amesos2

// namespace Ifpack2 {
// namespace Details {
//   // extern void __attribute__((weak)) registerLinearSolverFactory ();
//   // void __attribute__((weak)) registerLinearSolverFactory ();
//   // evoid __attribute__((weak)) registerLinearSolverFactory ();
// } // namespace Details
// } // namespace Ifpack2
#endif // ! defined(HAVE_TEUCHOS_DYNAMIC_LIBS) && defined(HAVE_TEUCHOS_CXX_ATTRIBUTE_WEAK)


/// \namespace Trilinos
/// \brief Namespace of things generally useful to many Trilinos packages
namespace Trilinos {

/// \namespace Details
/// \brief Namespace of implementation details.
///
/// \warning This namespace, and anything in it, is an implementation
///   detail of Trilinos.  Do not rely on this namespace or its
///   contents.  They may change or disappear at any time.
namespace Details {

template<class MV, class OP, class NormType>
class LinearSolver; // forward declaration

/// \brief Get a LinearSolver instance.
///
/// \tparam MV Type of a (multi)vector, representing either the
///   solution(s) X or the right-hand side(s) B of a linear system
///   AX=B.  For example, with Tpetra, use a Tpetra::MultiVector
///   specialization.  A <i>multivector</i> is a single data structure
///   containing zero or more vectors with the same dimensions and
///   layout.
///
/// \tparam OP Type of a matrix or linear operator that this Solver
///   understands.  For example, for Tpetra, use a Tpetra::Operator
///   specialization.
///
/// \tparam NormType Type of the norm of the residual.  See the
///   documentation of LinearSolver for details.
///
/// Call this function to create a LinearSolver instance from a
/// particular package.  LinearSolvers may create LinearSolvers.  The
/// run-time registration system (see registerLinearSolverFactory()
/// below) breaks software dependencies between packages.  Thus,
/// Package A may create a LinearSolver from Package B, even if
/// Package B depends on Package A.
///
/// \param packageName [in] Name of the package from which to get the
///   solver.  Names are case sensitive.
/// \param solverName [in] The solver's name.  Names are case sensitive.
template<class MV, class OP, class NormType>
Teuchos::RCP<LinearSolver<MV, OP, NormType> >
getLinearSolver (const std::string& packageName, const std::string& solverName);

/// \class LinearSolverFactory
/// \brief Interface for a "factory" that creates solvers.
///
/// \tparam MV Type of a (multi)vector, representing either the
///   solution(s) X or the right-hand side(s) B of a linear system
///   AX=B.  For example, with Tpetra, use a Tpetra::MultiVector
///   specialization.  A <i>multivector</i> is a single data structure
///   containing zero or more vectors with the same dimensions and
///   layout.
///
/// \tparam OP Type of a matrix or linear operator that the
///   LinearSolver to create understands.  For example, for Tpetra,
///   use a Tpetra::Operator specialization.
///
/// \tparam NormType Type of the norm of the residual.  See the
///   documentation of LinearSolver for details.
///
/// Every package that implements solvers needs to implement a
/// concrete LinearSolverFactory subclass as well.  That subclass
/// knows how to create all the solvers which that package implements.
/// The package must register a LinearSolverFactory instance using
/// registerLinearSolverFactory() (see below).  Then, any package may
/// access that package's solvers, using getLinearSolver() (see above).
///
/// You do not need to worry about "de-registering" or deallocating
/// LinearSolverFactory instances; std::shared_ptr takes care of that
/// automatically, after main() finishes.  LinearSolverFactory
/// instances should not hold on to resources that need explicit
/// deallocation or "finalization," such as MPI_* data structures
/// (that need to be "freed" before MPI_Finalize() is called) or open
/// file handles.
///
/// If you have a compelling use case that requires explicit
/// finalization of a LinearSolverFactory instance at some point
/// before main() finishes, please talk to the Trilinos developers
/// about adding a deregisterLinearSolverFactory() function (which
/// does not exist yet).
///
/// In the Tpetra solver stack, it's necessary to register factories
/// for all combinations of template parameters that applications plan
/// to use.  The easiest way to do that is to hook into the explicit
/// template instantiation (ETI) system of each package.  If ETI is
/// ON, this is easy.  If ETI is OFF, it's a bit harder.  Tpetra
/// defines a set of template parameter combinations over which it
/// _tests_.  If ETI is ON, this is always a subset of the ETI set.
/// If ETI is OFF, I would recommend using this set of test types for
/// registering factories.  Do the following:
///
/// 1. Include TpetraCore_ETIHelperMacros.h (a header file that
///    Tpetra's CMake configuration process generates and writes to
///    the build directory)
///
/// 2. In an anonymous outer namespace, define a class that registers
///    your factory in its constructor.  See PackageA.cpp in
///    ../example/SolverFactory for an example.  Then, define a macro
///    that creates an instance of that class (see the bottom of
///    PackageA.cpp).
///
/// 3. In the same anonymous outer namespace, invoke the
///    TPETRA_ETI_MANGLING_TYPEDEFS() macro to define typedefs used
///    internally by #3
///
/// 4. In the same anonymous outer namespace, use the
///    TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR macro, passing in the
///    name of your macro (see #2) as its one argument.
template<class MV, class OP, class NormType>
class LinearSolverFactory {
public:
  virtual ~LinearSolverFactory() {}
  /// \brief Get an instance of a solver from a particular package.
  ///
  /// \param solverName [in] The solver's name.  Names are case
  ///   sensitive.
  ///
  /// \return A pointer to the LinearSolver, if the name was valid;
  ///   else, a null pointer.
  virtual Teuchos::RCP<LinearSolver<MV, OP, NormType> >
  getLinearSolver (const std::string& solverName) = 0;
};

/// \function registerLinearSolverFactory
/// \brief Called by a package to register its LinearSolverFactory.
///
/// \note Most users do not need to call this function.  This is
///   mostly of interest to solver package developers.  See below for
///   details.
///
/// \tparam MV Type of a (multi)vector, representing either the
///   solution(s) X or the right-hand side(s) B of a linear system
///   AX=B.  For example, with Tpetra, use a Tpetra::MultiVector
///   specialization.  A <i>multivector</i> is a single data structure
///   containing zero or more vectors with the same dimensions and
///   layout.
///
/// \tparam OP Type of a matrix or linear operator that the
///   LinearSolver instances to create understand.  For example, for
///   Tpetra, use a Tpetra::Operator specialization.
///
/// \tparam NormType Type of the norm of the residual.  See the
///   documentation of LinearSolver for details.
///
/// \param packageName [in] Name of the package registering the
///   factory.  Package names are case sensitive.
/// \param factory [in] That package's factory.
///
/// This function lets packages register themselves, so that
/// getLinearSolver() (see above) can create solvers from that
/// package.  A package "registers itself" by doing the following:
/// <ol>
///   <li> Defining a concrete LinearSolverFactory subclass,
///        that knows how to create solvers from that package </li>
///   <li> Calling registerLinearSolverFactory() (this function) with
///        an instance of that LinearSolverFactory subclass </li>
/// </ol>
///
/// Packages may call this function before main() runs.  In fact, we
/// prefer that they do so.  This ensures that any package will be
/// able to create solvers from that package, without users or other
/// packages needing to know about that package.  When people talk
/// about "dependency injection" or "dependency inversion," this is
/// what they mean.
///
/// This function is templated with the same template parameters as
/// LinearSolverFactory.  This means that it must be called for every
/// combination of types (MV, OP) for which code will instantiate a
/// LinearSolverFactory<MV, OP, NormType>.  Thus, if the solver package wants to
/// do this before main() runs, it needs a list of all type
/// combination in advance.  If using explicit template instantiation
/// (ETI), you may plug this into the ETI system.  We thus recommend
/// that packages that use ETI register a LinearSolverFactory instance
/// for each ETI type combination.  For example, Ifpack2 should
/// iterate over all enabled combinations of the four template
/// parameters S, LO, GO, NT of Ifpack2::Preconditioner, creating a
/// LinearSolverFactory<MV, OP, NormType> instance for each combination, with MV
/// = Tpetra::MultiVector<S, LO, GO, NT> and OP = Tpetra::Operator<S,
/// LO, GO, NT>.  Package developers may find it useful to write a
/// macro that does this for that package's LinearSolverFactory
/// subclass.
///
/// If packages do not register a factory for certain type
/// combinations that users need, users may in rare instances need to
/// call this function themselves.  Avoid doing this, because it
/// defeats dependency inversion.
///
/// It could very well be that some packages don't implement all
/// desired type combinations MV, OP.  In that case, those packages
/// would not register a factory for those types.  Users who request
/// solvers from those packages for forbidden type combinations would
/// get a run-time error.
///
/// \note To developers: LinearSolverFactory returns LinearSolver by
///   Teuchos::RCP because Trilinos' solvers tend to use Teuchos::RCP,
///   and we don't want to break compatibility.  However, if C++11 is
///   enabled, we use std::shared_ptr to handle LinearSolverFactory
///   instances.  This is because that is an implementation detail
///   that solvers themselves don't have to see, and because
///   std::shared_ptr is thread safe.
template<class MV, class OP, class NormType>
void
registerLinearSolverFactory (const std::string& packageName,
#ifdef HAVE_TEUCHOSCORE_CXX11
                             const std::shared_ptr<LinearSolverFactory<MV, OP, NormType> >& factory);
#else
                             const Teuchos::RCP<LinearSolverFactory<MV, OP, NormType> >& factory);
#endif // HAVE_TEUCHOSCORE_CXX11

//
// EVERYTHING BELOW THIS LINE IS AN IMPLEMENTATION DETAIL
//

/// \brief Implementation details of implementation details.
///
/// We've already warned you that the Details namespace is full of
/// implementation details.  This inner namespace has implementation
/// details of <i>implementation details</i>.
namespace Impl {

/// \brief Remember which packages registered at least one
///   LinearSolverFactory, with any template parameters.
///
/// This is helpful for debugging failures to register a
/// LinearSolverFactory with the correct template parameters.
///
/// \return true if the package has already registered some
///   LinearSolverFactory before, else false.  (Same as what
///   registeredSomeLinearSolverFactory(packageName) would have
///   returned.)
bool rememberRegisteredSomeLinearSolverFactory (const std::string& packageName);

/// \brief Did the package with the given name register at least
///   one LinearSolverFactory, with any template parameters?
///
/// This is helpful for debugging failures to register a
/// LinearSolverFactory with the correct template parameters.
bool registeredSomeLinearSolverFactory (const std::string& packageName);

/// \brief Whether the CMake run-time registration option is ON.
///
/// This doesn't actually say whether run-time registration has
/// happened for a particular combination of (MV, OP, NormType)
/// template parameters.  Also, some packages or users may have
/// registered a factory manually; this has nothing to do with that.
bool haveLinearSolverFactoryRunTimeRegistration ();

/// \class LinearSolverFactoryRepository
/// \brief Repository of solver factories
///
/// \tparam MV Type of a (multi)vector, representing either the
///   solution(s) X or the right-hand side(s) B of a linear system
///   AX=B.  For example, with Tpetra, use a Tpetra::MultiVector
///   specialization.  A <i>multivector</i> is a single data structure
///   containing zero or more vectors with the same dimensions and
///   layout.
///
/// \tparam OP Type of a matrix or linear operator that LinearSolver
///   understands.  For example, for Tpetra, use a Tpetra::Operator
///   specialization.
///
/// \tparam NormType Type of the norm of a residual.
///
/// A LinearSolver knows how to solve linear systems AX=B.  A
/// LinearSolverFactory knows how to create LinearSolver instances.
/// Each independent unit of code ("package") that wants to
/// participate in the linear solver system, registers its own
/// LinearSolverFactory using the nonmember functions
/// Trilinos::Details::registerLinearSolverFactory().  Solvers may
/// then get (inner) solver instances with
/// Trilinos::Details::getLinearSolver() (see above in this file).
/// Those two nonmember functions dispatch to this class' class
/// (static) methods with the same names.
template<class MV, class OP, class NormType>
class LinearSolverFactoryRepository {
public:
  /// \typedef factory_pointer_type
  /// \brief Type of a reference-counted pointer to LinearSolverFactory.
  ///
  /// If C++11 is enabled, we use std::shared_ptr here, for improved
  /// thread safety.  Teuchos does not require C++11.
#ifdef HAVE_TEUCHOSCORE_CXX11
  typedef std::shared_ptr<LinearSolverFactory<MV, OP, NormType> > factory_pointer_type;
#else
  typedef Teuchos::RCP<LinearSolverFactory<MV, OP, NormType> > factory_pointer_type;
#endif // HAVE_TEUCHOSCORE_CXX11

  /// \typedef map_type
  /// \brief Type of a data structure that looks up a
  ///   LinearSolverFactory corresponding to a given package name.
  ///
  /// The compiler insists that this be public.  This doesn't hurt
  /// encapsulation, because this class lives in an "Impl"(ementation)
  /// namespace anyway.
  typedef std::map<std::string, factory_pointer_type> map_type;

public:
  /// \brief Get a LinearSolverFactory from the given package.
  ///
  /// This is an implementation detail of the nonmember function with
  /// the same name (see above).
  ///
  /// \param packageName [in] Name of the package.  This must be the
  ///   same name as that used to register the package via
  ///   registerLinearSolverFactory().  Package names are case
  ///   sensitive.
  ///
  /// \return If \c packageName has been registered with a valid
  ///   LinearSolverFactory, the pointer to the factory, else null.
  static factory_pointer_type
  getFactory (const std::string& packageName)
  {
    createFactories ();
    typedef typename map_type::iterator iter_type;
    iter_type it = factories_->find (packageName);
    if (it == factories_->end ()) { // didn't find package name
      return factory_pointer_type (); // null pointer
    } else { // found package name
      return it->second;
    }
  }

  /// \brief Register the given factory from a package.
  ///
  /// This is an implementation detail of the nonmember function with
  /// the same name (see above).
  ///
  /// \param packageName [in] Name of the package registering the
  ///   factory.  Package names are case sensitive.
  /// \param factory [in] That package's factory (must be nonnull).
  ///
  /// \warning This method is not reentrant.  In particular, if
  ///   multiple threads call this method at the same time, they might
  ///   manage to double-register the atexit() handler for factories_.
  ///   This could only happen if this method is called twice by
  ///   different threads.
  static void
  registerLinearSolverFactory (const std::string& packageName,
                               const factory_pointer_type& factory)
  {
    TEUCHOS_TEST_FOR_EXCEPTION
      (factory.get () == NULL, std::invalid_argument, "Trilinos::Details::"
       "LinearSolverFactoryRepository::registerLinearSolverFactory: Input "
       "'factory' is NULL!");
    createFactories ();
    if (factories_->find (packageName) == factories_->end ()) {
      factories_->insert (std::make_pair (packageName, factory));
    }
  }

private:
  /// \brief Singleton where all packages' factories get stored.
  ///
  /// The map maps from each package's name (as given to
  /// registerLinearSolverFactory()) to a pointer to that package's
  /// LinearSolverFactory instance.
  ///
  /// This unfortunately has to be a pointer.  Otherwise, the std::map
  /// never gets initialized, and segfaults result.  We initialize the
  /// pointer in createFactories(), where we set an atexit() hook to
  /// free it using freeFactories().
  static map_type* factories_;

  /// \brief Initialize factories_ if it hasn't been initialized.
  ///
  /// Also, set an atexit() hook to free it using freeFactories().
  static void createFactories () {
    if (factories_ == NULL) {
      factories_ = new map_type ();
      // It _is_ possible for atexit() to fail (e.g., because it ran
      // out of memory for storing callbacks).  We could throw an
      // exception here in that case, but I think it's better just
      // to let the minor memory leak happen.
      (void) atexit (freeFactories);
    }
    TEUCHOS_TEST_FOR_EXCEPTION
      (factories_ == NULL, std::logic_error, "Trilinos::Details::"
       "LinearSolverFactoryRepository::createFactories: "
       "Should never get here!  factories_ is NULL.");
  }

  /// \brief Free the factories_ singleton.
  ///
  /// \warning Only for use as atexit() handler.
  ///
  /// \warning This method is not reentrant.  In particular, if
  ///   multiple threads call this method at the same time, they might
  ///   manage to double-delete factories_.  This should not happen
  ///   because the atexit() hook should only ever be called once.
  static void freeFactories () {
    if (factories_ != NULL) {
      delete factories_;
      factories_ = NULL;
    }
  }
};

// This is _not_ an explicit instantiation.  C++ wants it, because
// LinearSolverFactoryRepository is a templated class with a static
// (class) member.
template<class MV, class OP, class NormType>
typename LinearSolverFactoryRepository<MV, OP, NormType>::map_type*
LinearSolverFactoryRepository<MV, OP, NormType>::factories_ = NULL;

} // namespace Impl

//
// Definitions of nonmember functions
//

template<class MV, class OP, class NormType>
void
registerLinearSolverFactory (const std::string& packageName,
#ifdef HAVE_TEUCHOSCORE_CXX11
                             const std::shared_ptr<LinearSolverFactory<MV, OP, NormType> >& factory)
#else
                             const Teuchos::RCP<LinearSolverFactory<MV, OP, NormType> >& factory)
#endif // HAVE_TEUCHOSCORE_CXX11
{
  Impl::LinearSolverFactoryRepository<MV, OP, NormType>::registerLinearSolverFactory (packageName, factory);
  Impl::rememberRegisteredSomeLinearSolverFactory (packageName);
}

template<class MV, class OP, class NormType>
Teuchos::RCP<LinearSolver<MV, OP, NormType> >
getLinearSolver (const std::string& packageName, const std::string& solverName)
{
  using Teuchos::RCP;
  using Teuchos::TypeNameTraits;
  typedef Impl::LinearSolverFactoryRepository<MV, OP, NormType> repo_type;
  typedef typename repo_type::factory_pointer_type factory_pointer_type;
  typedef LinearSolver<MV, OP, NormType> solver_type;
  const char prefix[] = "Trilinos::Details::getLinearSolver: ";

  // FIXME (mfh 21 Aug 2015) Attempted fix for Bug 6392: DOES NOT WORK.
  // (Compiles just fine, but test doesn't pass.)
#if ! defined(HAVE_TEUCHOS_DYNAMIC_LIBS) && defined(HAVE_TEUCHOS_CXX_ATTRIBUTE_WEAK)
  // if (Amesos2::Details::registerLinearSolverFactory == NULL) {
  //   std::cout << "-- Amesos2::Details::registerLinearSolverFactory is NULL" << std::endl;
  // } else {
  //   Amesos2::Details::registerLinearSolverFactory ();
  // }
  // if (Ifpack2::Details::registerLinearSolverFactory == NULL) {
  //   std::cout << "-- Ifpack2::Details::registerLinearSolverFactory is NULL" << std::endl;
  // } else {
  //   Ifpack2::Details::registerLinearSolverFactory ();
  // }
#endif // ! defined(HAVE_TEUCHOS_DYNAMIC_LIBS) && defined(HAVE_TEUCHOS_CXX_ATTRIBUTE_WEAK)

  // Whether the CMake run-time registration option is ON.  This
  // doesn't actually say whether run-time registration has happened
  // for the current combination of (MV, OP, NormType) template
  // parameters.
  const bool haveRunTimeReg =
    Impl::haveLinearSolverFactoryRunTimeRegistration ();

  const bool pkgExists = Impl::registeredSomeLinearSolverFactory (packageName);
  TEUCHOS_TEST_FOR_EXCEPTION
    (! pkgExists, std::invalid_argument, prefix << "Package \"" << packageName
     << "\" never registered a LinearSolverFactory for _any_ combination of "
     "template parameters MV, OP, and NormType.  This means either that the "
     "package name is invalid, or that the package is not enabled.  "
     "Trilinos_ENABLE_LINEAR_SOLVER_FACTORY_REGISTRATION = "
     << (haveRunTimeReg ? "ON" : "OFF") << ".");

  factory_pointer_type factory = repo_type::getFactory (packageName);
  TEUCHOS_TEST_FOR_EXCEPTION
    (factory.get () == NULL, std::invalid_argument, prefix << "Package \"" <<
     packageName << "\" is valid, but it never registered a LinearSolverFactory"
     " for template parameters "
     "MV = " << TypeNameTraits<MV>::name () << ", "
     "OP = " << TypeNameTraits<OP>::name () << ", "
     "NormType = " << TypeNameTraits<NormType>::name () << ".  "
     "Trilinos_ENABLE_LINEAR_SOLVER_FACTORY_REGISTRATION = "
     << (haveRunTimeReg ? "ON" : "OFF") << ".");

  RCP<solver_type> solver = factory->getLinearSolver (solverName);
  TEUCHOS_TEST_FOR_EXCEPTION
    (solver.is_null (), std::invalid_argument, prefix << "Invalid solver name "
     "\"" << solverName << "\".  However, package \"" << packageName << "\" is "
     "valid, and it did register a LinearSolverFactory for template parameters "
     "MV = " << TypeNameTraits<MV>::name () << ", "
     "OP = " << TypeNameTraits<OP>::name () << ", "
     "NormType = " << TypeNameTraits<NormType>::name () << ".  "
     "Trilinos_ENABLE_LINEAR_SOLVER_FACTORY_REGISTRATION = "
     << (haveRunTimeReg ? "ON" : "OFF") << ".");

  return solver;
}

} // namespace Details
} // namespace Trilinos

#endif // TRILINOS_DETAILS_LINEARSOLVER_FACTORY_HPP
