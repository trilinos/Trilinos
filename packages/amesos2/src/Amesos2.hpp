// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Amesos2_config.h"
#include "Amesos2_Factory.hpp"

// Define some groups to be used in our doxygen documentation

/**
 * \defgroup amesos2_adapters Amesos2 Linear Algebra Object Adapters
 *
 * Amesos2 derives much of its flexibility from its adapters.  Amesos2
 * solver instances are templated on matrix and multivector types.  As
 * long as an Amesos2 adapter exists for those types, then Amesos2 can
 * interact with those objects.
 *
 * Amesos2 has two types of adapters:
 * - \ref amesos2_matrix_adapters "Matrix Adapters", and
 * - \ref amesos2_multivec_adapters "MultiVector Adapters"
 *
 * The adapters provide a unifying interface for Amesos2 to use,
 * regardless of the actual type of the object.  In this way, a solver
 * interface itself does not need to be changed in order to work with
 * a new linear algebra object.
 * 
 * New adapters can be created whenever the need arises.  Just contact
 * the Amesos2 developers if you have a linear algebra object that you
 * want Amesos2 to be able to interact with.
 */

/**
 * \defgroup amesos2_matrix_adapters Amesos2 Matrix Adapters
 * \ingroup amesos2_adapters
 *
 * Amesos2 matrix adapters provide an interface that caters to many of
 * the existing third-party sparse direct solver libraries.  Most (if
 * not all) third-party libraries accept matrices in either the
 * compressed sparse row format or the compressed sparse column
 * format.  Retrieving such a representation from a matrix is the
 * major task of the Amesos2::MatrixAdapter class (along with
 * providing other critical statistics about a matrix, such as
 * row/column size, number of nonzero entries, etc).
 *
 * The Amesos2::AbstractConcreteMatrixAdapter and
 * Amesos2::ConcreteMatrixAdapter, as a pair, exist to fully exploit
 * inheritence in the adaptation.  The ConcreteMatrixAdapter is where
 * the MatrixAdapter class goes to for much of its business.  For
 * example, when getting a compressed sparse row representation of a
 * matrix, MatrixAdapter will depend on the ConcreteMatrixAdapter's
 * getGlobalRowCopy() method.  However, for both the Epetra and Tpetra
 * matrices, there typically exists some sort of abstraction, such as
 * the Epetra_RowMatrix interface.  The Amesos matrix adapters are
 * able to exploit such interfaces (abstractions) through the
 * AbstractConcreteMatrixAdapter template class.  The purpose of each
 * specialization of this templated class is to adapt as many of the
 * abstract methods (those methods specified in the abstract
 * interface) as possible.  The ConcreteMatrixAdapter can then inherit
 * much of the adapted functionality as possible.  Take for instance
 * the Epetra_RowMatrix abstraction.  This abstraction specifies a set
 * of methods that deriving classes must implement (or simply inherit
 * the default implementations).  A number of Epetra matrices derive
 * from this abstraction.  By adapting the abstraction, rather than
 * the implementing classes, the Amesos2 adapters significantly reduce
 * redundant code.
 */

/**
 * \defgroup amesos2_multivec_adapters Amesos2 MultiVector Adapters
 * \ingroup amesos2_adapters
 *
 * Amesos2 solver interfaces are most interested in being able to
 * access and update the values found in multivectors.  This is the
 * primary goal of the Amesos2 multivector adapters.  They provide
 * methods for querying properties of the multivector (number of
 * vectors, length, etc), for getting copies of the values stored in
 * the multivector, and for updating the values (e.g. once a solution
 * vector has been found, put this solution in the multivector).
 */

/////////////////////////// Solvers ////////////////////////////////////

/**
 * \defgroup amesos2_solvers Amesos2 Solvers
 *
 * Perhaps the most interesting part of Amesos2 from a user's
 * perspective, but the Amesos2 solver system was designed to be
 * useful for both users and developers.  The system can be split into
 * two distinct but inter-related parts: The \ref
 * amesos2_solver_framework "solver framework", and the \ref
 * amesos2_solver_interfaces "solver interfaces".
 *
 * The Amesos2::Solver class provides a uniform interface to the
 * third-party library solvers.  The interface is designed to be both
 * simple to use for novice users, as well as powerful enough for
 * advanced users.  While a novice user might like to just give a
 * linear system to Amesos2 and have it just solve it, an expert user
 * might like to control how and when each step of the solution
 * process is performed and do solves for multiple different RHS
 * vectors.
 *
 * An example of solving a system with Amesos2 using it's most basic
 * interface:
 *
 * \code
 * RCP<MAT> A; RCP<MV> X; RCP<MV> B;
 * // initialize A and B
 * RCP<Solver<MAT,MV> > solver = Amesos2::create(A, X, B); // use default solver
 * solver->solve(); // solution placed in X
 * \endcode
 *
 * Here is another more involved example:
 *
 * \code
 * RCP<MAT> A;
 * // Get A from somewhere
 * RCP<Solver<MAT,MV> > solver = Amesos2::create("SuperLU", A);
 * Teuchos::ParameterList params("Amesos2");
 * params.sublist("SuperLU").set("IterRefine","DOUBLE");
 * params.sublist("SuperLU").set("ColPerm","MMD_AT_PLUS_A");
 * solver->setParameters(params);
 * solver->symbolicFactorization().numericFactorization();
 * A = Teuchos::null;          // no longer need A
 * solver.setA(Teuchos::null); // tell the solver to release A too
 * RCP<MV> X; RCP<MV> B;
 * // do some other work, finally get B's values
 * solver->solve(X,B);	       // solution placed in X
 * // do some more work and get new values in B
 * solver->solve(X,B);
 * \endcode
 */

/**
 * \defgroup amesos2_solver_framework Amesos2 Solver Framework
 * \ingroup amesos2_solvers
 *
 * The Amesos2 Solver Framework provides an infrastructure that \ref
 * amesos2_solver_interfaces "concrete solver interfaces" depend upon.
 * At the same time this framework provides a sort of
 * fill-in-the-blank framework for solving a system of linear
 * equations that depends upon the solver-interfaces to fill in those
 * blanks.
 * 
 * Its purpose is to abstract all the solver-independent features in
 * order to reduce the burden on those developing new solver
 * interfaces.  In doing this it reduces the amount of maintained code
 * and focuses a developers concerns.  A developer writing a new
 * solver interface does not have to be bothered to recreate
 * infrastructure but can instead rely on the framework provided by
 * Amesos2.
 *
 * We could describe the framework in terms of the "Chain of
 * Responsibility" pattern.  When a user requests for a solve to be
 * done, the SolveCore class takes care of any work that would need to
 * be done by any solver (such as starting some timers, checking
 * for exceptions, updating internal state, etc.) and then passes of
 * responsibility to the concrete solver interfaces, who in turn
 * delegate the algorithmic work to our third-party libraries.
 *
 * \note In a way, the \ref amesos2_adapters "Amesos2 adapters" could
 * be viewed as part of the solver framework, but they are interesting
 * enough in their own right to be listed as a separate module.
 */

/**
 * \defgroup amesos2_solver_interfaces Amesos2 Solver Interfaces
 * \ingroup amesos2_solvers
 *
 * The Amesos2 solver interfaces are responsible for distilling a
 * third-party library's interface into an Amesos2-like interface.
 * For the most part, the solver interfaces need only concern
 * themselves with getting data to the TPL and storing results for
 * later or sending solutions back to the user.
 *
 * \note 
 * Users of Amesos2 do not need to concern themselves directly with
 * the solver interface classes
 */

//////////////////// Solver Parameters ////////////////////

/**
 * \defgroup amesos2_solver_parameters Supported Solver Parameters
 * \ingroup amesos2_solvers
 *
 * Many third-party solvers support a vast amount of parameters to
 * control factorization and solution.  An effort has been made in
 * Amesos2 to expose to users as many of those parameters as
 * reasonably possible.  Not all parameters may be supported, but if
 * there is one that you would like to have exposed, then contact the
 * Amesos2 developers and we may be able to work something out for
 * you.
 *
 * \section amesos2_parameters
 *
 * The following parameters are currently acted upon by Amesos2
 * solvers:
 * 
 * <ul>
 *   <li> \c "Transpose" : { \c true | \c false }.  If \c true , tells
 *   the solver to solve for \f$A^T X = B\f$</li>
 * </ul>
 *
 * We plan in the future to support the following parameters:
 *
 * <ul>
 *   <li> \c "Reindex" : { \c true | \c false }.  Put the matrix row
 *     and column indices into the range [0..n].</li>
 *   <li> \c "AddZeroToDiag" : { \c true | \c false }.</li>
 *   <li> \c "AddToDiag" : \c string . Where the given string is a
 *     representation of a scalar value to add to all diagonal
 *     entries of the matrix before factorization.</li>
 * </ul>
 *
 * \section amesos2_solver_parameters Solver-specific Parameters
 *
 * \subsection superlu_parameters SuperLU
 *
 * \copydoc Amesos2::Superlu::setParameters_impl()
 *
 * \subsection superlu_mt_parameters SuperLU_MT
 *
 * \copydetails Amesos2::Superlumt::setParameters_impl()
 *
 * \subsection superlu_dist_parameters SuperLU_DIST
 *
 * \copydetails Amesos2::Superludist::setParameters_impl()
 *
 * \subsection pardiso_mkl_parameters Pardiso MKL
 *
 * \copydetails Amesos2::PardisoMKL::setParameters_impl()
 *
 * \subsection lapack_parameters LAPACK
 *
 * \copydetails Amesos2::Lapack::setParameters_impl()
 *
 * \subsection umfpack_parameters Umfpack
 *
 * \copydetails Amesos2::Umfpack::setParameters_impl()
 */


//////////////////// Types ////////////////////

/**
 * \defgroup amesos2_enums Amesos2 Enum Types
 */


//////////////////////////////// Examples /////////////////////////////////

/*
 * We put all the example blocks in a single place, since doxygen doesn't do
 * an incredible job with putting links to the examples where you would like
 * them.
 */

/**
 * \example SimpleSolve.cpp
 *
 * Shows how to create an Amesos2 solver using the Amesos2::create()
 * factory method interface, followed by solving a small linear
 * system.
 */

/**
 * \example SimpleSolve_File.cpp
 *
 * Shows how you could use Amesos2 with Tpetra's MatrixMarket
 * functionality.
 */

/**
 * \example SimpleSolve_WithParameters.cpp
 *
 * An example of how to set solver parameters with the Superlu interface.
 */

/**
 * \example MultipleSolves_File.cpp
 *
 * This example gives an example of how to re-solve a system after
 * having changed it, either changing the matrix itself, or changing
 * part of the RHS B.
 */

/**
 * \example TwoPartSolve.cpp
 *
 * An example of how one can defer providing the X and B vectors to
 * the Amesos2 solver until just before calling a solve using the setX
 * and setB methods.
 */
