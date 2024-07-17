// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file   Amesos2_Details_LinearSolverFactory_def.hpp
/// \author Mark Hoemmen
/// \brief  Definition of Amesos2::Details::LinearSolverFactory.

#ifndef AMESOS2_DETAILS_LINEARSOLVERFACTORY_DEF_HPP
#define AMESOS2_DETAILS_LINEARSOLVERFACTORY_DEF_HPP

#include "Amesos2_config.h"
#include "Amesos2_Factory.hpp"
#include "Amesos2_Solver.hpp"
#include "Trilinos_Details_LinearSolver.hpp"
#include "Trilinos_Details_LinearSolverFactory.hpp"
#include <type_traits>


#ifdef HAVE_AMESOS2_EPETRA
#  include "Epetra_CrsMatrix.h"
#endif // HAVE_AMESOS2_EPETRA

// mfh 23 Jul 2015: Tpetra is currently a required dependency of Amesos2.
#ifndef HAVE_AMESOS2_TPETRA
#  define HAVE_AMESOS2_TPETRA
#endif // HAVE_AMESOS2_TPETRA

#ifdef HAVE_AMESOS2_TPETRA
#  include "Tpetra_CrsMatrix.hpp"
#endif // HAVE_AMESOS2_TPETRA

namespace Amesos2 {
namespace Details {

// For a given linear algebra implementation's Operator type OP,
// find the corresponding CrsMatrix type.
//
// Amesos2 unfortunately only does ETI for Tpetra::CrsMatrix, even
// though it could very well take Tpetra::RowMatrix.
template<class OP>
struct GetMatrixType {
  typedef int type; // flag (see below)


#ifdef HAVE_AMESOS2_EPETRA
  static_assert(! std::is_same<OP, Epetra_MultiVector>::value,
                "Amesos2::Details::GetMatrixType: OP = Epetra_MultiVector.  "
                "This probably means that you mixed up MV and OP.");
#endif // HAVE_AMESOS2_EPETRA

#ifdef HAVE_AMESOS2_TPETRA
  static_assert(! std::is_same<OP, Tpetra::MultiVector<typename OP::scalar_type,
                typename OP::local_ordinal_type, typename OP::global_ordinal_type,
                typename OP::node_type> >::value,
                "Amesos2::Details::GetMatrixType: OP = Tpetra::MultiVector.  "
                "This probably means that you mixed up MV and OP.");
#endif // HAVE_AMESOS2_TPETRA
};

#ifdef HAVE_AMESOS2_EPETRA
template<>
struct GetMatrixType<Epetra_Operator> {
  typedef Epetra_CrsMatrix type;
};
#endif // HAVE_AMESOS2_EPETRA


#ifdef HAVE_AMESOS2_TPETRA
template<class S, class LO, class GO, class NT>
struct GetMatrixType<Tpetra::Operator<S, LO, GO, NT> > {
  typedef Tpetra::CrsMatrix<S, LO, GO, NT> type;
};
#endif // HAVE_AMESOS2_TPETRA

template<class MV, class OP, class NormType>
class LinearSolver :
    public Trilinos::Details::LinearSolver<MV, OP, NormType>,
    virtual public Teuchos::Describable
{

#ifdef HAVE_AMESOS2_EPETRA
  static_assert(! std::is_same<OP, Epetra_MultiVector>::value,
                "Amesos2::Details::LinearSolver: OP = Epetra_MultiVector.  "
                "This probably means that you mixed up MV and OP.");
  static_assert(! std::is_same<MV, Epetra_Operator>::value,
                "Amesos2::Details::LinearSolver: MV = Epetra_Operator.  "
                "This probably means that you mixed up MV and OP.");
#endif // HAVE_AMESOS2_EPETRA

public:
  /// \brief Type of the CrsMatrix specialization corresponding to OP.
  ///
  /// Amesos2 currently requires a CrsMatrix as its template
  /// parameter, not an Operator or a RowMatrix.  However, the
  /// SolverFactory system is templated on Operator.  This type traits
  /// class converts from OP to the corresponding CrsMatrix type.
  typedef typename GetMatrixType<OP>::type crs_matrix_type;
  static_assert(! std::is_same<crs_matrix_type, int>::value,
                "Amesos2::Details::LinearSolver: The given OP type is not "
                "supported.");

  //! Type of the underlying Amesos2 solver.
  typedef Amesos2::Solver<crs_matrix_type, MV> solver_type;

  /// \brief Constructor.
  ///
  /// \param solverName [in] Name of the solver to create.  Must be a
  ///   name accepted by Amesos2::create.
  ///
  /// \note To developers: The constructor takes a name rather than a
  ///   solver, because Amesos2 does not currently let us create a
  ///   solver with a null matrix.  To get around this, we create a
  ///   new solver whenever the matrix changes.  The unfortunate side
  ///   effect is that, since Amesos2 has no way of determining
  ///   whether a solver supports a given set of template parameters
  ///   without attempting to create the solver, we can't actually
  ///   test in the constructor whether the solver exists.
  LinearSolver (const std::string& solverName) :
    solverName_ (solverName)
  {
    // FIXME (mfh 25 Aug 2015) Ifpack2::Details::Amesos2Wrapper made
    // the unfortunate choice to attempt to guess solvers that exist.
    // Thus, alas, for the sake of backwards compatibility, we must
    // make an attempt to guess for the user.  Furthermore, for strict
    // backwards compatibility, we should preserve the same (rather
    // arbitrary) list of choices, in the same order.
    if (solverName == "") {
      // KLU2 is the default solver.
      if (Amesos2::query ("klu2")) {
        solverName_ = "klu2";
      }
      else if (Amesos2::query ("superlu")) {
        solverName_ = "superlu";
      }
      else if (Amesos2::query ("superludist")) {
        solverName_ = "superludist";
      }
      else if (Amesos2::query ("cholmod")) {
        solverName_ = "cholmod";
      }
      else if (Amesos2::query ("cusolver")) {
        solverName_ = "cusolver";
      }
      else if (Amesos2::query ("basker")) {
        solverName_ = "basker";
      }
      else if (Amesos2::query ("shylubasker")) {
        solverName_ = "shylubasker";
      }
      else if (Amesos2::query ("ShyLUBasker")) {
        solverName_ = "shylubasker";
      }
      else if (Amesos2::query ("superlumt")) {
        solverName_ = "superlumt";
      }
      else if (Amesos2::query ("pardiso_mkl")) {
        solverName_ = "pardiso_mkl";
      }
      else if (Amesos2::query ("css_mkl")) {
        solverName_ = "css_mkl";
      }
      else if (Amesos2::query ("mumps")) {
        solverName_ = "mumps";
      }
      else if (Amesos2::query ("lapack")) {
        solverName_ = "lapack";
      }
      else if (Amesos2::query ("umfpack")) {
        solverName_ = "umfpack";
      }
      // We don't have to try to rescue the user if their empty solver
      // name didn't catch any of the above choices.
    }
  }

  //! Destructor (virtual for memory safety).
  virtual ~LinearSolver () {}

  /// \brief Set the Solver's matrix.
  ///
  /// \param A [in] Pointer to the matrix A in the linear system(s)
  ///   AX=B to solve.
  void setMatrix (const Teuchos::RCP<const OP>& A) {
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::TypeNameTraits;
    typedef crs_matrix_type MAT;
    const char prefix[] = "Amesos2::Details::LinearSolver::setMatrix: ";

    if (A.is_null ()) {
      solver_ = Teuchos::null;
    }
    else {
      // FIXME (mfh 24 Jul 2015) Even though Amesos2 solvers currently
      // require a CrsMatrix input, we could add a copy step in order
      // to let them take a RowMatrix.  The issue is that this would
      // require keeping the original matrix around, and copying again
      // if it changes.
      RCP<const MAT> A_mat = Teuchos::rcp_dynamic_cast<const MAT> (A);
      TEUCHOS_TEST_FOR_EXCEPTION
        (A_mat.is_null (), std::invalid_argument,
         "Amesos2::Details::LinearSolver::setMatrix: "
         "The input matrix A must be a CrsMatrix.");
      if (solver_.is_null ()) {
        // Amesos2 solvers must be created with a nonnull matrix.
        // Thus, we don't actually create the solver until setMatrix
        // is called for the first time with a nonnull matrix.
        // Furthermore, Amesos2 solvers do not accept a null matrix
        // (to setA), so if setMatrix is called with a null matrix, we
        // free the solver.
        RCP<solver_type> solver;
        try {
          solver = Amesos2::create<MAT, MV> (solverName_, A_mat, null, null);
        }
        catch (std::exception& e) {
          TEUCHOS_TEST_FOR_EXCEPTION
            (true, std::invalid_argument, prefix << "Failed to create Amesos2 "
             "solver named \"" << solverName_ << "\".  "
             "Amesos2::create<MAT = " << TypeNameTraits<MAT>::name ()
             << ", MV = " << TypeNameTraits<MV>::name ()
             << " threw an exception: " << e.what ());
        }
        TEUCHOS_TEST_FOR_EXCEPTION
          (solver.is_null (), std::invalid_argument, prefix << "Failed to "
           "create Amesos2 solver named \"" << solverName_ << "\".  "
           "Amesos2::create<MAT = " << TypeNameTraits<MAT>::name ()
           << ", MV = " << TypeNameTraits<MV>::name ()
           << " returned null.");

        // Use same parameters as before, if user set parameters.
        if (! params_.is_null ()) {
          solver->setParameters (params_);
        }
        solver_ = solver;
      } else if (A_ != A) {
        solver_->setA (A_mat);
      }
    }

    // We also have to keep a pointer to A, so that getMatrix() works.
    A_ = A;
  }

  //! Get a pointer to this Solver's matrix.
  Teuchos::RCP<const OP> getMatrix () const {
    return A_;
  }

  //! Solve the linear system(s) AX=B.
  void solve (MV& X, const MV& B) {
    const char prefix[] = "Amesos2::Details::LinearSolver::solve: ";
    TEUCHOS_TEST_FOR_EXCEPTION
      (solver_.is_null (), std::runtime_error, prefix << "The solver does not "
       "exist yet.  You must call setMatrix() with a nonnull matrix before you "
       "may call this method.");
    TEUCHOS_TEST_FOR_EXCEPTION
      (A_.is_null (), std::runtime_error, prefix << "The matrix has not been "
       "set yet.  You must call setMatrix() with a nonnull matrix before you "
       "may call this method.");
    solver_->solve (&X, &B);
  }

  //! Set this solver's parameters.
  void setParameters (const Teuchos::RCP<Teuchos::ParameterList>& params) {
    if (! solver_.is_null ()) {
      solver_->setParameters (params);
    }
    // Remember them, so that if the solver gets recreated, we'll have
    // the original parameters.
    params_ = params;
  }

  /// \brief Set up any part of the solve that depends on the
  ///   structure of the input matrix, but not its numerical values.
  void symbolic () {
    const char prefix[] = "Amesos2::Details::LinearSolver::symbolic: ";
    TEUCHOS_TEST_FOR_EXCEPTION
      (solver_.is_null (), std::runtime_error, prefix << "The solver does not "
       "exist yet.  You must call setMatrix() with a nonnull matrix before you "
       "may call this method.");
    TEUCHOS_TEST_FOR_EXCEPTION
      (A_.is_null (), std::runtime_error, prefix << "The matrix has not been "
       "set yet.  You must call setMatrix() with a nonnull matrix before you "
       "may call this method.");
    solver_->symbolicFactorization ();
  }

  /// \brief Set up any part of the solve that depends on both the
  ///   structure and the numerical values of the input matrix.
  void numeric () {
    const char prefix[] = "Amesos2::Details::LinearSolver::numeric: ";
    TEUCHOS_TEST_FOR_EXCEPTION
      (solver_.is_null (), std::runtime_error, prefix << "The solver does not "
       "exist yet.  You must call setMatrix() with a nonnull matrix before you "
       "may call this method.");
    TEUCHOS_TEST_FOR_EXCEPTION
      (A_.is_null (), std::runtime_error, prefix << "The matrix has not been "
       "set yet.  You must call setMatrix() with a nonnull matrix before you "
       "may call this method.");
    solver_->numericFactorization ();
  }

  //! Implementation of Teuchos::Describable::description.
  std::string description () const {
    using Teuchos::TypeNameTraits;
    if (solver_.is_null ()) {
      std::ostringstream os;
      os << "\"Amesos2::Details::LinearSolver\": {"
         << "MV: " << TypeNameTraits<MV>::name ()
         << ", OP: " << TypeNameTraits<OP>::name ()
         << ", NormType: " << TypeNameTraits<NormType>::name ()
         << "}";
      return os.str ();
    } else {
      return solver_->description ();
    }
  }

  //! Implementation of Teuchos::Describable::describe.
  void
  describe (Teuchos::FancyOStream& out,
            const Teuchos::EVerbosityLevel verbLevel =
            Teuchos::Describable::verbLevel_default) const
  {
    using Teuchos::TypeNameTraits;
    using std::endl;
    if (solver_.is_null ()) {
      if (verbLevel > Teuchos::VERB_NONE) {
        Teuchos::OSTab tab0 (out);
        out << "\"Amesos2::Details::LinearSolver\":" << endl;
        Teuchos::OSTab tab1 (out);
        out << "MV: " << TypeNameTraits<MV>::name () << endl
            << "OP: " << TypeNameTraits<OP>::name () << endl
            << "NormType: " << TypeNameTraits<NormType>::name () << endl;
      }
    }
    if (! solver_.is_null ()) {
      solver_->describe (out, verbLevel);
    }
  }

private:
  std::string solverName_;
  Teuchos::RCP<solver_type> solver_;
  Teuchos::RCP<const OP> A_;
  Teuchos::RCP<Teuchos::ParameterList> params_;
};

template<class MV, class OP, class NormType>
Teuchos::RCP<Trilinos::Details::LinearSolver<MV, OP, NormType> >
LinearSolverFactory<MV, OP, NormType>::
getLinearSolver (const std::string& solverName)
{
  using Teuchos::rcp;
  return rcp (new Amesos2::Details::LinearSolver<MV, OP, NormType> (solverName));
}

template<class MV, class OP, class NormType>
void
LinearSolverFactory<MV, OP, NormType>::
registerLinearSolverFactory ()
{
#ifdef HAVE_TEUCHOSCORE_CXX11
  typedef std::shared_ptr<Amesos2::Details::LinearSolverFactory<MV, OP, NormType> > ptr_type;
  //typedef std::shared_ptr<Trilinos::Details::LinearSolverFactory<MV, OP> > base_ptr_type;
#else
  typedef Teuchos::RCP<Amesos2::Details::LinearSolverFactory<MV, OP, NormType> > ptr_type;
  //typedef Teuchos::RCP<Trilinos::Details::LinearSolverFactory<MV, OP> > base_ptr_type;
#endif // HAVE_TEUCHOSCORE_CXX11

  ptr_type factory (new Amesos2::Details::LinearSolverFactory<MV, OP, NormType> ());
  Trilinos::Details::registerLinearSolverFactory<MV, OP, NormType> ("Amesos2", factory);
}

} // namespace Details
} // namespace Amesos2

// Macro for doing explicit instantiation of
// Amesos2::Details::LinearSolverFactory, for Tpetra objects, with
// given Tpetra template parameters (SC = Scalar, LO = LocalOrdinal,
// GO = GlobalOrdinal, NT = Node).
//
// We don't have to protect use of Tpetra objects here, or include
// any header files for them, because this is a macro definition.
#define AMESOS2_DETAILS_LINEARSOLVERFACTORY_INSTANT(SC, LO, GO, NT) \
  template class Amesos2::Details::LinearSolverFactory<Tpetra::MultiVector<SC, LO, GO, NT>, \
                                                       Tpetra::Operator<SC, LO, GO, NT>, \
                                                       typename Tpetra::MultiVector<SC, LO, GO, NT>::mag_type>;

#endif // AMESOS2_DETAILS_LINEARSOLVERFACTORY_DEF_HPP
