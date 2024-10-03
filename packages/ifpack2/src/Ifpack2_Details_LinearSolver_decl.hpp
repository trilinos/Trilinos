// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Ifpack2_Details_LinearSolver_decl_hpp
/// \brief Declaration of Ifpack2::Details::LinearSolver, an
///   implementation detail of Ifpack2's LinearSolverFactory.

#ifndef IFPACK2_DETAILS_LINEARSOLVER_DECL_HPP
#define IFPACK2_DETAILS_LINEARSOLVER_DECL_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Trilinos_Details_LinearSolver.hpp"
#include "Ifpack2_Preconditioner.hpp"
#include "Teuchos_Describable.hpp"

namespace Ifpack2 {
namespace Details {

/// \class LinearSolver
/// \brief Ifpack2's implementation of Trilinos::Details::LinearSolver
///   interface
///
/// \tparam SC Scalar type; 1st template parameter of Tpetra::Operator
/// \tparam LO Local ordinal type; 2nd template parameter of Tpetra::Operator
/// \tparam GO Global ordinal type; 3rd template parameter of Tpetra::Operator
/// \tparam NT Node type; 4th template parameter of Tpetra::Operator
///
/// Note to Ifpack2 developers: The main choice when implementing
/// Trilinos::Details::LinearSolver for a solver package, is whether
/// the LinearSolver or the LinearSolverFactory should know how to
/// create solvers from that package.  For Ifpack2, one would think
/// that LinearSolver needs to know how to create solvers.  This is
/// because not every Ifpack2::Preconditioner subclass knows how to
/// change its matrix, which is necessary for implementing setMatrix.
/// If solvers don't know how to change their matrix, then the
/// LinearSolver subclass has to destroy and recreate the solver with
/// the new matix.
///
/// Ifpack2 has a complication, though: Some Ifpack2 solvers have a
/// choice, either to create an inner solver using
/// Trilinos::Details::LinearSolverFactory, or to use an
/// Ifpack2::Preconditioner provided by the user.  AdditiveSchwarz is
/// the prototypical example.  For simplicity, we let AdditiveSchwarz
/// wrap the user's Ifpack2::Preconditioner in this LinearSolver
/// class.  However, that means that AdditiveSchwarz needs to include
/// Ifpack2_Details_LinearSolver.hpp (since it needs to invoke the
/// constructor).  When ETI (explicit template instantiation) is OFF,
/// this means that AdditiveSchwarz pulls in both the declaration and
/// definition of this LinearSolver class.  As a result, we cannnot
/// use Ifpack2::Factory in this class, because that would reintroduce
/// one of the circular dependencies that LinearSolverFactory is meant
/// to avoid!
///
/// Thankfully, Ifpack2 gives us a way out.  AdditiveSchwarz already
/// requires that its inner solver implement
/// Ifpack2::Details::CanChangeMatrix, which has a setMatrix() method.
/// This means we can restrict the set of Ifpack2 solvers that
/// LinearSolver wraps, to those that implement CanChangeMatrix.  As a
/// result, we can move use of Ifpack2::Factory to
/// Ifpack2::Details::LinearSolverFactory.  This class thus takes an
/// Ifpack2::Preconditioner, and a solver name.  If the solver name is
/// "CUSTOM", we assume that the Ifpack2::Preconditioner comes from
/// the user (via e.g., AdditiveSchwarz::setInnerPreconditioner), and
/// that Ifpack2::Factory might not know how to create it.
template<class SC, class LO, class GO, class NT>
class LinearSolver :
    public Trilinos::Details::LinearSolver<Tpetra::MultiVector<SC, LO, GO, NT>,
                                           Tpetra::Operator<SC, LO, GO, NT>,
                                           typename Tpetra::MultiVector<SC, LO, GO, NT>::mag_type>,
    virtual public Teuchos::Describable
{
public:
  typedef Ifpack2::Preconditioner<SC, LO, GO, NT> prec_type;
  typedef Tpetra::Operator<SC, LO, GO, NT> OP;
  typedef Tpetra::MultiVector<SC, LO, GO, NT> MV;

  /// \brief Constructor
  ///
  /// \param solver [in] The Ifpack2 solver to wrap.  It MUST
  ///   implement Ifpack2::Details::CanChangeMatrix<row_matrix_type>.
  ///
  /// \param solverName [in] Name of the solver.  If the solver name
  ///   is "CUSTOM", we assume that the solver comes from the user
  ///   (via e.g., Ifpack2::AdditiveSchwarz::setInnerPreconditioner),
  ///   and that Ifpack2::Factory might not know how to create it.
  ///   Otherwise, the name needs to be that of a solver that
  ///   Ifpack2::Factory::create knows how to create.
  LinearSolver (const Teuchos::RCP<prec_type>& solver, const std::string& solverName);

  //! Destructor (virtual for memory safety).
  virtual ~LinearSolver () {}

  /// \brief Set the solver's matrix.
  ///
  /// \param A [in] Pointer to the matrix A in the linear system(s)
  ///   AX=B to solve.
  void setMatrix (const Teuchos::RCP<const OP>& A);

  //! Get the solver's matrix.
  Teuchos::RCP<const OP> getMatrix () const;

  //! Solve the linear system AX=B for X.
  void solve (MV& X, const MV& B);

  //! Set the solver's parameters.
  void setParameters (const Teuchos::RCP<Teuchos::ParameterList>& params);

  //! Precompute for matrix structure changes.
  void symbolic ();

  //! Precompute for matrix values' changes.
  void numeric ();

  //! Implementation of Teuchos::Describable::description.
  std::string description () const;

  //! Implementation of Teuchos::Describable::describe.
  void
  describe (Teuchos::FancyOStream& out,
            const Teuchos::EVerbosityLevel verbLevel =
            Teuchos::Describable::verbLevel_default) const;

private:
  //! The Ifpack2 solver.
  Teuchos::RCP<prec_type> solver_;
  //! Name of the Ifpack2 solver to wrap.
  std::string solverName_;
  //! Matrix A in the linear system to solve.
  Teuchos::RCP<const OP> A_;
};

} // namespace Details
} // namespace Ifpack2

#endif // IFPACK2_DETAILS_LINEARSOLVER_DECL_HPP
