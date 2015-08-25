#ifndef IFPACK2_DETAILS_LINEARSOLVER_DECL_HPP
#define IFPACK2_DETAILS_LINEARSOLVER_DECL_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Trilinos_Details_LinearSolver.hpp"
#include "Ifpack2_Preconditioner.hpp"
#include "Teuchos_Describable.hpp"

namespace Ifpack2 {
namespace Details {

/// \class LinearSolver
/// \brief Ifpack2's implementation of Trilinos::Details::LinearSolver interface
///
/// \tparam SC Scalar type; 1st template parameter of Tpetra::Operator
/// \tparam LO Local ordinal type; 2nd template parameter of Tpetra::Operator
/// \tparam GO Global ordinal type; 3rd template parameter of Tpetra::Operator
/// \tparam NT Node type; 4th template parameter of Tpetra::Operator
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
  typedef Tpetra::RowMatrix<SC, LO, GO, NT> row_matrix_type;
  typedef Tpetra::MultiVector<SC, LO, GO, NT> MV;

  /// \brief Constructor
  ///
  /// \param solverName [in] Name of the solver to create.  Must be a
  ///   name accepted by Ifpack2::Factory::create.
  ///
  /// \note To developers: The constructor takes a name rather than a
  ///   solver, because Ifpack2 does not currently let us create a
  ///   solver with a null matrix.  To get around this, we create a
  ///   new solver whenever the matrix changes.  The unfortunate side
  ///   effect is that, since Ifpack2 has no way of determining
  ///   whether a solver supports a given set of template parameters
  ///   without attempting to create the solver, we can't actually
  ///   test in the constructor whether the solver exists.
  LinearSolver (const std::string& solverName);

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
  //! Name of the Ifpack2 solver to wrap.
  std::string solverName_;
  //! The Ifpack2 solver.
  Teuchos::RCP<prec_type> solver_;
  //! Matrix A in the linear system to solve.
  Teuchos::RCP<const OP> A_;
  //! The solver's parameters.
  Teuchos::RCP<Teuchos::ParameterList> params_;
};

} // namespace Details
} // namespace Ifpack2

#endif // IFPACK2_DETAILS_LINEARSOLVER_DECL_HPP
