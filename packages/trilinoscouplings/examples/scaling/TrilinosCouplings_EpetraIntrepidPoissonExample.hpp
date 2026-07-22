// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __TrilinosCouplings_EpetraIntrepidPoissonExample_hpp
#define __TrilinosCouplings_EpetraIntrepidPoissonExample_hpp

#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_ScalarTraits.hpp"


namespace TrilinosCouplings {
/// \namespace EpetraIntrepidPoissonExample
/// \brief Epetra version of the Intrepid Poisson test problem example.
///
/// The Intrepid Poisson test problem uses Pamgen to construct a 3-D
/// mesh (a simple rectangular prism with hex elements) in parallel,
/// Sacado automatic differentiation to construct a right-hand side of
/// the PDE corresponding to a given exact solution, and Intrepid to
/// build a discretization.
///
/// We provide two variants of the Intrepid Poisson test: one that
/// fills Epetra objects, and one that fills Tpetra objects.  The two
/// variants do exactly the same things otherwise, so you can use them
/// to compare the performance of Epetra and Tpetra fill.
///
/// This namespace contains the Epetra variant.  It defines typedefs
/// which you can use when writing a main() driver to run the test.
/// The makeMatrixAndRightHandSide() function does all the work.  You
/// can use the exactResidualNorm() function to test correctness of
/// the discretization.  In particular, if the continuous exact
/// solution was chosen from the space of finite element polynomials,
/// the exact solution of the discrete linear system AX=B should match
/// the continuous exact solution exactly, modulo rounding error when
/// assembling the discretization.
///
/// The solveWithBelos() function solves the given linear system using
/// a Belos iterative solver.  You can provide a left and/or right
/// preconditioner if you want.
namespace EpetraIntrepidPoissonExample {

//
// Typedefs for use by main() and other functions.
//
typedef double ST;
typedef int    LO;
typedef int    GO;
typedef Epetra_CrsMatrix   sparse_matrix_type;
typedef Epetra_Operator    operator_type;
typedef Epetra_MultiVector multivector_type;
typedef Epetra_Vector      vector_type;

/// \brief Create the mesh and build the linear system to solve.
///
/// \param A [out] The sparse matrix.
/// \param B [out] The right-hand side(s).
/// \param X_exact [out] The exact solution of the PDE, projected onto
///   the discrete mesh.  This may not necessarily be the same as the
///   exact solution of the discrete linear system.
/// \param X [out] The approximate solution(s).
/// \param err [out] Output stream for errors.
/// \param out [out] Output stream for verbose output.
/// \param comm [in] Communicator (as an Epetra_Comm).
/// \param meshInput [in] Pamgen mesh specification string.
///
/// See the documentation of
/// TpetraIntrepidPoissonExample::makeMatrixAndRightHandSide() for an
/// explanation of Pamgen and a reference for the mesh specification
/// format.
void
makeMatrixAndRightHandSide (Teuchos::RCP<sparse_matrix_type>& A,
                            Teuchos::RCP<vector_type>& B,
                            Teuchos::RCP<vector_type>& X_exact,
                            Teuchos::RCP<vector_type>& X,
                            const Teuchos::RCP<const Epetra_Comm>& comm,
                            const std::string& meshInput,
                            const Teuchos::RCP<Teuchos::FancyOStream>& out,
                            const Teuchos::RCP<Teuchos::FancyOStream>& err,
                            const bool verbose = false,
                            const bool debug = false);

//! Just like above, but with multivector_type output arguments.
void
makeMatrixAndRightHandSide (Teuchos::RCP<sparse_matrix_type>& A,
                            Teuchos::RCP<multivector_type>& B,
                            Teuchos::RCP<multivector_type>& X_exact,
                            Teuchos::RCP<multivector_type>& X,
                            const Teuchos::RCP<const Epetra_Comm>& comm,
                            const std::string& meshInput,
                            const Teuchos::RCP<Teuchos::FancyOStream>& out,
                            const Teuchos::RCP<Teuchos::FancyOStream>& err,
                            const bool verbose = false,
                            const bool debug = false);

//! Return \f$\|B - A X_{\text{exact}}\|_2\f$, \f$\|B\|\f$, and \f$\|A\|_F\f$.
std::vector<Teuchos::ScalarTraits<ST>::magnitudeType>
exactResidualNorm (const Teuchos::RCP<const sparse_matrix_type>& A,
                   const Teuchos::RCP<const vector_type>& B,
                   const Teuchos::RCP<const vector_type>& X_exact);

/// \brief Solve the linear system(s) AX=B with Belos.
///
/// X and B are both multivectors, meaning that you may ask Belos to
/// solve more than one linear system at a time.  X and B must have
/// the same number of columns.
///
/// This interface will change in the future to accept the name of the
/// Belos solver to use.  For now, the solver is hard-coded to
/// Pseudoblock CG (implemented by Belos::PseudoBlockCGSolMgr).
///
/// \param converged [out] Whether Belos reported that the iterative
///   method converged (solved the linear system to the desired
///   tolerance).
///
/// \param numItersPerformed [out] Number of iterations that the Belos
///   solver performed.
///
/// \param solverName [in] Which iterative linear solver to use.
///   Any name that Belos::SolverFactory knows will work here.
///
/// \param tol [in] Convergence tolerance for the iterative method.
///   The meaning of this depends on the particular iterative method.
///
/// \param maxNumIters [in] Maximum number of iterations that the
///   iterative method should perform, regardless of whether it
///   converged.
///
/// \param num_steps [in] Number of "time steps", i.e., the number of
//    times the solver is called in a fake time-step loop.
///
/// \param X [in/out] On input: the initial guess(es) for the iterative
///   method.  On output: the computed approximate solution.
///
/// \param A [in] The matrix in the linear system(s) AX=B to solve.
///
/// \param B [in] The right-hand side(s) in the linear system AX=B to solve.
///
/// \param M_left [in] If nonnull, a left preconditioner that the
///   iterative method may use.  If null, the iterative method will
///   not use a left preconditioner.
///
/// \param M_right [in] If nonnull, a right preconditioner that the
///   iterative method may use.  If null, the iterative method will
///   not use a right preconditioner.
void
solveWithBelos (bool& converged,
                int& numItersPerformed,
                const std::string& solverName,
                const Teuchos::ScalarTraits<ST>::magnitudeType& tol,
                const int maxNumIters,
                const int num_steps,
                const Teuchos::RCP<multivector_type>& X,
                const Teuchos::RCP<const sparse_matrix_type>& A,
                const Teuchos::RCP<const multivector_type>& B,
                const Teuchos::RCP<operator_type>& M_left=Teuchos::null,
                const Teuchos::RCP<operator_type>& M_right=Teuchos::null);

} // namespace EpetraIntrepidPoissonExample
} // namespace TrilinosCouplings

#endif // __TrilinosCouplings_EpetraIntrepidPoissonExample_hpp
