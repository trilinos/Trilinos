// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOSSOLVERFACTORYBASE_HPP
#define BELOSSOLVERFACTORYBASE_HPP

#include <string>
#include <vector>
#include "Teuchos_RCP.hpp"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Teuchos {
// Forward declaration
class ParameterList;
} // namespace Teuchos
#endif // DOXYGEN_SHOULD_SKIP_THIS

namespace Belos {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
// Forward declaration
template<class Scalar, class MV, class OP>
class SolverManager;
#endif // DOXYGEN_SHOULD_SKIP_THIS

/// \brief Interface for custom Belos solver factories
///
/// If you want to extend Belos::SolverFactory (which see) to support
/// new solvers, you may do the following:
///
/// <ol>
/// <li> Create a subclass of CustomSolverFactory, whose getSolver()
///      method returns instances of the new solvers </li>
/// <li> Tell Belos about your subclass, by giving an instance of it
///      to Belos::SolverFactory::addFactory </li>
/// <li> Now, you may use Belos::SolverFactory to create instances of
///      the new solvers </li>
/// </ol>
///
/// For a test and example of how to do this (with a trivial solver),
/// see <tt>Trilinos/packages/belos/tpetra/test/CustomSolverFactory.cpp</tt>.
template<class Scalar, class MV, class OP>
class CustomSolverFactory {
public:
  /// \brief Return an instance of the specified solver, or
  ///   Teuchos::null if this factory does not provide the requested
  ///   solver.
  ///
  /// \note To implementers: DO NOT THROW if this factory does not
  ///   recognize and support the input solverName.  Instead, just
  ///   return Teuchos::null.  This will make Belos::SolverFactory's
  ///   implementation cleaner (and is also the reason why we gave
  ///   this function a different name).
  ///
  /// \param solverName [in] Name of the requested solver.
  ///
  /// \param solverParams [in/out] List of parameters with which to
  ///   configure the solver.  If Teuchos::null, subclasses must still
  ///   return a valid solver.  They may configure the solver with
  ///   default parameters in that case.  If nonnull, the factory or
  ///   solver may modify the list by filling in missing parameters,
  ///   or overriding existing parameters.  You may then inspect the
  ///   resulting list to learn what parameters the solver accepts.
  ///
  /// \return If the given solverName is valid, the solver, else
  ///   Teuchos::null.  (Returning null makes "chaining" factories
  ///   work; see e.g., Belos::SolverFactory::addFactory.)
  ///
  /// The input parameter list is passed in as a Teuchos::RCP because
  /// the factory passes it to the solver, and Belos solvers want
  /// their input parameter list as a
  /// Teuchos::RCP<Teuchos::ParameterList>.  We allow a null parameter
  /// list only for convenience, and will use default parameter values
  /// in that case.
  virtual Teuchos::RCP<SolverManager<Scalar, MV, OP> >
  getSolver (const std::string& solverName,
             const Teuchos::RCP<Teuchos::ParameterList>& solverParams) = 0;

  /// \brief Number of supported solvers.
  ///
  /// This may differ from the number of supported solver
  /// <i>names</i>, since we may accept multiple names ("aliases") for
  /// some solvers.
  virtual int numSupportedSolvers () const = 0;

  /// \brief List of supported solver names.
  ///
  /// The length of this list may differ from the number of supported
  /// solvers, since we may accept multiple names ("aliases") for some
  /// solvers.
  virtual std::vector<std::string> supportedSolverNames () const = 0;

  //! Whether the given solver name names a supported solver.
  virtual bool isSupported (const std::string& solverName) const = 0;

  //! Destructor (virtual, for safety of derived classes).
  virtual ~CustomSolverFactory () {}
};

} // namespace Belos

#endif // BELOSSOLVERFACTORYBASE_HPP
