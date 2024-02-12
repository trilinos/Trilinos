// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2011 Sandia Corporation
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
// ***********************************************************************
//
// @HEADER

/// \file   MueLu_Details_LinearSolverFactory_def.hpp
/// \authors Mark Hoemmen and Alicia Klinvex
/// \brief  Definition of MueLu::Details::LinearSolverFactory.

#ifndef MUELU_DETAILS_LINEARSOLVERFACTORY_DEF_HPP
#define MUELU_DETAILS_LINEARSOLVERFACTORY_DEF_HPP

#include "MueLu_config.hpp"
#include "Trilinos_Details_LinearSolver.hpp"
#include "Trilinos_Details_LinearSolverFactory.hpp"
#include <type_traits>

#ifdef HAVE_MUELU_EPETRA
#include "Epetra_CrsMatrix.h"
#include "MueLu_CreateEpetraPreconditioner.hpp"
#endif  // HAVE_MUELU_EPETRA

#include "Tpetra_Operator.hpp"
#include "MueLu_CreateTpetraPreconditioner.hpp"

namespace MueLu {
namespace Details {

template <class MV, class OP, class NormType>
class LinearSolver : public Trilinos::Details::LinearSolver<MV, OP, NormType>,
                     virtual public Teuchos::Describable {
 public:
  /// \brief Constructor.
  LinearSolver() {}

  //! Destructor (virtual for memory safety).
  virtual ~LinearSolver() {}

  /// \brief Set the Solver's matrix.
  ///
  /// \param A [in] Pointer to the matrix A in the linear system(s)
  ///   AX=B to solve.
  void setMatrix(const Teuchos::RCP<const OP>& A);

  //! Get a pointer to this Solver's matrix.
  Teuchos::RCP<const OP> getMatrix() const {
    return A_;
  }

  //! Solve the linear system(s) AX=B.
  void solve(MV& X, const MV& B);

  //! Set this solver's parameters.
  void setParameters(const Teuchos::RCP<Teuchos::ParameterList>& params);

  /// \brief Set up any part of the solve that depends on the
  ///   structure of the input matrix, but not its numerical values.
  void symbolic() {}

  /// \brief Set up any part of the solve that depends on both the
  ///   structure and the numerical values of the input matrix.
  void numeric();

  //! Implementation of Teuchos::Describable::description.
  std::string description() const;

  //! Implementation of Teuchos::Describable::describe.
  void
  describe(Teuchos::FancyOStream& out,
           const Teuchos::EVerbosityLevel verbLevel =
               Teuchos::Describable::verbLevel_default) const;

 private:
  Teuchos::RCP<const OP> A_;
  Teuchos::RCP<Teuchos::ParameterList> params_;
};

// Why does MueLu_EpetraOperator insist on HAVE_MUELU_SERIAL?
#if defined(HAVE_MUELU_SERIAL) and defined(HAVE_MUELU_EPETRA)
template <>
class LinearSolver<Epetra_MultiVector, Epetra_Operator, double> : public Trilinos::Details::LinearSolver<Epetra_MultiVector, Epetra_Operator, double>,
                                                                  virtual public Teuchos::Describable {
 public:
  /// \brief Constructor.
  LinearSolver()
    : changedA_(false)
    , changedParams_(false) {}

  //! Destructor (virtual for memory safety).
  virtual ~LinearSolver() {}

  /// \brief Set the Solver's matrix.
  ///
  /// \param A [in] Pointer to the matrix A in the linear system(s)
  ///   AX=B to solve.
  void setMatrix(const Teuchos::RCP<const Epetra_Operator>& A) {
    const char prefix[] = "MueLu::Details::LinearSolver::setMatrix: ";

    if (A != A_) {
      if (solver_ != Teuchos::null)
        changedA_ = true;

      A_ = rcp_dynamic_cast<const Epetra_CrsMatrix>(A);
      TEUCHOS_TEST_FOR_EXCEPTION(A_.is_null(), std::runtime_error, prefix << "MueLu requires "
                                                                             "an Epetra_CrsMatrix, but the matrix you provided is of a "
                                                                             "different type.  Please provide an Epetra_CrsMatrix instead.");
    }
  }

  //! Get a pointer to this Solver's matrix.
  Teuchos::RCP<const Epetra_Operator> getMatrix() const {
    return A_;
  }

  //! Solve the linear system(s) AX=B.
  void solve(Epetra_MultiVector& X, const Epetra_MultiVector& B) {
    // TODO amk: Do we assume the user has called numeric before solve, or should we call it for them?
    const char prefix[] = "MueLu::Details::LinearSolver::solve: ";
    TEUCHOS_TEST_FOR_EXCEPTION(solver_.is_null(), std::runtime_error, prefix << "The solver does not "
                                                                                "exist yet.  You must call numeric() before you may call this method.");
    TEUCHOS_TEST_FOR_EXCEPTION(changedA_, std::runtime_error, prefix << "The matrix A has been reset "
                                                                        "since the last call to numeric().  Please call numeric() again.");
    TEUCHOS_TEST_FOR_EXCEPTION(changedParams_, std::runtime_error, prefix << "The parameters have been reset "
                                                                             "since the last call to numeric().  Please call numeric() again.");

    int err = solver_->ApplyInverse(B, X);

    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, prefix << "EpetraOperator::ApplyInverse returned "
                                                                       "nonzero error code "
                                                                    << err);
  }

  //! Set this solver's parameters.
  void setParameters(const Teuchos::RCP<Teuchos::ParameterList>& params) {
    if (solver_ != Teuchos::null && params != params_)
      changedParams_ = true;

    params_ = params;
  }

  /// \brief Set up any part of the solve that depends on the
  ///   structure of the input matrix, but not its numerical values.
  void symbolic() {}

  /// \brief Set up any part of the solve that depends on both the
  ///   structure and the numerical values of the input matrix.
  void numeric() {
    const char prefix[] = "MueLu::Details::LinearSolver::numeric: ";

    // If the solver is up-to-date, leave it alone
    if (solver_ == Teuchos::null || changedA_ || changedParams_) {
      changedA_      = false;
      changedParams_ = false;

      TEUCHOS_TEST_FOR_EXCEPTION(A_ == Teuchos::null, std::runtime_error, prefix << "The matrix has not been "
                                                                                    "set yet.  You must call setMatrix() with a nonnull matrix before you may "
                                                                                    "call this method.");

      // TODO: We should not have to cast away the constness here
      // TODO: See bug 6462
      if (params_ != Teuchos::null)
        solver_ = CreateEpetraPreconditioner(rcp_const_cast<Epetra_CrsMatrix>(A_), *params_);
      else
        solver_ = CreateEpetraPreconditioner(rcp_const_cast<Epetra_CrsMatrix>(A_));
    }
  }

  //! Implementation of Teuchos::Describable::description.
  std::string description() const {
    if (solver_.is_null()) {
      return "\"MueLu::Details::LinearSolver\": {MV: Epetra_MultiVector, OP: Epetra_Operator, NormType: double}";
    } else {
      return solver_->GetHierarchy()->description();
    }
  }

  //! Implementation of Teuchos::Describable::describe.
  void
  describe(Teuchos::FancyOStream& out,
           const Teuchos::EVerbosityLevel verbLevel =
               Teuchos::Describable::verbLevel_default) const {
    using std::endl;
    if (solver_.is_null()) {
      if (verbLevel > Teuchos::VERB_NONE) {
        Teuchos::OSTab tab0(out);
        out << "\"MueLu::Details::LinearSolver\":" << endl;
        Teuchos::OSTab tab1(out);
        out << "MV: Epetra_MultiVector" << endl
            << "OP: Epetra_Operator" << endl
            << "NormType: double" << endl;
      }
    } else {
      solver_->GetHierarchy()->describe(out, verbLevel);
    }
  }

 private:
  Teuchos::RCP<const Epetra_CrsMatrix> A_;
  Teuchos::RCP<Teuchos::ParameterList> params_;
  Teuchos::RCP<EpetraOperator> solver_;
  bool changedA_;
  bool changedParams_;
};
#endif  // HAVE_MUELU_EPETRA

template <class Scalar, class LO, class GO, class Node>
class LinearSolver<Tpetra::MultiVector<Scalar, LO, GO, Node>,
                   Tpetra::Operator<Scalar, LO, GO, Node>,
                   typename Teuchos::ScalarTraits<Scalar>::magnitudeType> : public Trilinos::Details::LinearSolver<Tpetra::MultiVector<Scalar, LO, GO, Node>,
                                                                                                                   Tpetra::Operator<Scalar, LO, GO, Node>,
                                                                                                                   typename Teuchos::ScalarTraits<Scalar>::magnitudeType>,
                                                                            virtual public Teuchos::Describable {
 public:
  /// \brief Constructor.
  LinearSolver()
    : changedA_(false)
    , changedParams_(false) {}

  //! Destructor (virtual for memory safety).
  virtual ~LinearSolver() {}

  /// \brief Set the Solver's matrix.
  ///
  /// \param A [in] Pointer to the matrix A in the linear system(s)
  ///   AX=B to solve.
  void setMatrix(const Teuchos::RCP<const Tpetra::Operator<Scalar, LO, GO, Node> >& A) {
    if (A != A_) {
      if (solver_ != Teuchos::null)
        changedA_ = true;

      A_ = A;
    }
  }

  //! Get a pointer to this Solver's matrix.
  Teuchos::RCP<const Tpetra::Operator<Scalar, LO, GO, Node> > getMatrix() const {
    return A_;
  }

  //! Solve the linear system(s) AX=B.
  void solve(Tpetra::MultiVector<Scalar, LO, GO, Node>& X, const Tpetra::MultiVector<Scalar, LO, GO, Node>& B) {
    // TODO amk: Do we assume the user has called numeric before solve, or should we call it for them?
    const char prefix[] = "MueLu::Details::LinearSolver::solve: ";
    TEUCHOS_TEST_FOR_EXCEPTION(solver_.is_null(), std::runtime_error, prefix << "The solver does not "
                                                                                "exist yet.  You must call numeric() before you may call this method.");
    TEUCHOS_TEST_FOR_EXCEPTION(changedA_, std::runtime_error, prefix << "The matrix A has been reset "
                                                                        "since the last call to numeric().  Please call numeric() again.");
    TEUCHOS_TEST_FOR_EXCEPTION(changedParams_, std::runtime_error, prefix << "The parameters have been reset "
                                                                             "since the last call to numeric().  Please call numeric() again.");

    solver_->apply(B, X);
  }

  //! Set this solver's parameters.
  void setParameters(const Teuchos::RCP<Teuchos::ParameterList>& params) {
    if (solver_ != Teuchos::null && params != params_)
      changedParams_ = true;

    params_ = params;
  }

  /// \brief Set up any part of the solve that depends on the
  ///   structure of the input matrix, but not its numerical values.
  void symbolic() {}

  /// \brief Set up any part of the solve that depends on both the
  ///   structure and the numerical values of the input matrix.
  void numeric() {
    const char prefix[] = "MueLu::Details::LinearSolver::numeric: ";

    // If the solver is up-to-date, leave it alone
    if (solver_ == Teuchos::null || changedParams_) {
      TEUCHOS_TEST_FOR_EXCEPTION(A_ == Teuchos::null, std::runtime_error, prefix << "The matrix has not been "
                                                                                    "set yet.  You must call setMatrix() with a nonnull matrix before you may "
                                                                                    "call this method.");

      // TODO: We should not have to cast away the constness here
      // TODO: See bug 6462
      if (params_ != Teuchos::null)
        solver_ = CreateTpetraPreconditioner(rcp_const_cast<Tpetra::Operator<Scalar, LO, GO, Node> >(A_), *params_);
      else
        solver_ = CreateTpetraPreconditioner(rcp_const_cast<Tpetra::Operator<Scalar, LO, GO, Node> >(A_));
    } else if (changedA_) {
      TEUCHOS_TEST_FOR_EXCEPTION(A_ == Teuchos::null, std::runtime_error, prefix << "The matrix has not been "
                                                                                    "set yet.  You must call setMatrix() with a nonnull matrix before you may "
                                                                                    "call this method.");

      // TODO: We should not have to cast away the constness here
      // TODO: See bug 6462
      RCP<const Tpetra::CrsMatrix<Scalar, LO, GO, Node> > helperMat;
      helperMat = rcp_dynamic_cast<const Tpetra::CrsMatrix<Scalar, LO, GO, Node> >(A_);
      TEUCHOS_TEST_FOR_EXCEPTION(helperMat.is_null(), std::runtime_error, prefix << "MueLu requires "
                                                                                    "a Tpetra::CrsMatrix, but the matrix you provided is of a "
                                                                                    "different type.  Please provide a Tpetra::CrsMatrix instead.");
      ReuseTpetraPreconditioner(rcp_const_cast<Tpetra::CrsMatrix<Scalar, LO, GO, Node> >(helperMat), *solver_);
    }

    changedA_      = false;
    changedParams_ = false;
  }

  //! Implementation of Teuchos::Describable::description.
  std::string description() const {
    using Teuchos::TypeNameTraits;
    if (solver_.is_null()) {
      std::ostringstream os;
      os << "\"MueLu::Details::LinearSolver\": {"
         << "MV: " << TypeNameTraits<Tpetra::MultiVector<Scalar, LO, GO, Node> >::name()
         << "OP: " << TypeNameTraits<Tpetra::Operator<Scalar, LO, GO, Node> >::name()
         << "NormType: " << TypeNameTraits<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>::name()
         << "}";
      return os.str();
    } else {
      return solver_->GetHierarchy()->description();
    }
  }

  //! Implementation of Teuchos::Describable::describe.
  void
  describe(Teuchos::FancyOStream& out,
           const Teuchos::EVerbosityLevel verbLevel =
               Teuchos::Describable::verbLevel_default) const {
    using std::endl;
    using Teuchos::TypeNameTraits;
    if (solver_.is_null()) {
      if (verbLevel > Teuchos::VERB_NONE) {
        Teuchos::OSTab tab0(out);
        out << "\"MueLu::Details::LinearSolver\":" << endl;
        Teuchos::OSTab tab1(out);
        out << "MV: " << TypeNameTraits<Tpetra::MultiVector<Scalar, LO, GO, Node> >::name() << endl
            << "OP: " << TypeNameTraits<Tpetra::Operator<Scalar, LO, GO, Node> >::name() << endl
            << "NormType: " << TypeNameTraits<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>::name() << endl;
      }
    } else {
      solver_->GetHierarchy()->describe(out, verbLevel);
    }
  }

 private:
  Teuchos::RCP<const Tpetra::Operator<Scalar, LO, GO, Node> > A_;
  Teuchos::RCP<Teuchos::ParameterList> params_;
  Teuchos::RCP<TpetraOperator<Scalar, LO, GO, Node> > solver_;
  bool changedA_;
  bool changedParams_;
};

template <class MV, class OP, class NormType>
Teuchos::RCP<Trilinos::Details::LinearSolver<MV, OP, NormType> >
LinearSolverFactory<MV, OP, NormType>::
    getLinearSolver(const std::string& solverName) {
  using Teuchos::rcp;
  return rcp(new MueLu::Details::LinearSolver<MV, OP, NormType>());
}

template <class MV, class OP, class NormType>
void LinearSolverFactory<MV, OP, NormType>::
    registerLinearSolverFactory() {
#ifdef HAVE_TEUCHOSCORE_CXX11
  typedef std::shared_ptr<MueLu::Details::LinearSolverFactory<MV, OP, NormType> > ptr_type;
  // typedef std::shared_ptr<Trilinos::Details::LinearSolverFactory<MV, OP> > base_ptr_type;
#else
  typedef Teuchos::RCP<MueLu::Details::LinearSolverFactory<MV, OP, NormType> > ptr_type;
  // typedef Teuchos::RCP<Trilinos::Details::LinearSolverFactory<MV, OP> > base_ptr_type;
#endif  // HAVE_TEUCHOSCORE_CXX11

  ptr_type factory(new MueLu::Details::LinearSolverFactory<MV, OP, NormType>());
  Trilinos::Details::registerLinearSolverFactory<MV, OP, NormType>("MueLu", factory);
}

}  // namespace Details
}  // namespace MueLu

// Macro for doing explicit instantiation of
// MueLu::Details::LinearSolverFactory, for Tpetra objects, with
// given Tpetra template parameters (SC = Scalar, LO = LocalOrdinal,
// GO = GlobalOrdinal, NT = Node).
//
// We don't have to protect use of Tpetra objects here, or include
// any header files for them, because this is a macro definition.
#define MUELU_DETAILS_LINEARSOLVERFACTORY_INSTANT(SC, LO, GO, NT)                         \
  template class MueLu::Details::LinearSolverFactory<Tpetra::MultiVector<SC, LO, GO, NT>, \
                                                     Tpetra::Operator<SC, LO, GO, NT>,    \
                                                     typename Tpetra::MultiVector<SC, LO, GO, NT>::mag_type>;

#endif  // MUELU_DETAILS_LINEARSOLVERFACTORY_DEF_HPP
