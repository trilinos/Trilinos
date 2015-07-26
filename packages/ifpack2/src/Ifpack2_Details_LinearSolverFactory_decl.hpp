/// \file   Ifpack2_Details_LinearSolverFactory_decl.hpp
/// \author Mark Hoemmen
/// \brief  Declaration of Ifpack2::Details::LinearSolverFactory.

#ifndef IFPACK2_DETAILS_LINEARSOLVERFACTORY_DECL_HPP
#define IFPACK2_DETAILS_LINEARSOLVERFACTORY_DECL_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Trilinos_Details_LinearSolverFactory.hpp"

namespace Ifpack2 {
namespace Details {

  /// \class LinearSolverFactory
  /// \brief Interface for a "factory" that creates Ifpack2 solvers.
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
  ///   specialization.  Always use the most abstract interface
  ///   possible; solvers should dynamic_cast to the subclass they
  ///   need.  Also, be consistent: using different classes here
  ///   (e.g., Tpetra::RowMatrix instead of Tpetra::Operator) means
  ///   more expensive explicit template instantiation.
  template<class MV, class OP>
  class LinearSolverFactory :
    public Trilinos::Details::LinearSolverFactory<MV, OP> {
  public:
    /// \brief Get an instance of a Ifpack2 solver.
    ///
    /// The solver is wrapped in a Trilinos::Details::LinearSolver
    /// interface.
    ///
    /// \param solverName [in] The solver's name.  Names are case
    ///   sensitive.
    /// \return A pointer to the solver, if the name was valid; else,
    ///   a null pointer (Teuchos::null).
    virtual Teuchos::RCP<Trilinos::Details::LinearSolver<MV, OP> >
    getLinearSolver (const std::string& solverName);
  };

} // namespace Details
} // namespace Ifpack2

#endif // IFPACK2_DETAILS_LINEARSOLVERFACTORY_DECL_HPP
