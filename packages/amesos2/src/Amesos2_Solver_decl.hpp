// @HEADER
// ***********************************************************************
//
//                Amesos2: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

/**
  \file   Amesos2_Solver_decl.hpp
  \author Eric T Bavier <etbavier@sandia.gov>
  \date   Thu May 27 13:51:05 CDT 2010

  \brief  Templated class for Amesos2 solvers
*/

#ifndef AMESOS2_SOLVER_DECL_HPP
#define AMESOS2_SOLVER_DECL_HPP

#include <string>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "Amesos2_SolverBase.hpp"
#include "Amesos2_MatrixAdapter.hpp"
#include "Amesos2_MultiVecAdapter.hpp"
#include "Amesos2_Util.hpp"


namespace Amesos {


/* This is the base class to be used in a *statically* polymorphic
 * way. E.g. for the Superlu solver:
 *
 * In Amesos2_Superlu.hpp:
 * class Superlu : Solver<Superlu> { ... }
 *
 * Each concrete solver will implement several private subfunction
 * that will be called within the common code for each function.
 *
 * In addition, the Amesos2::Solver static base class provides a convenient
 * way to handle the template parameters and the private typedef-ing (with
 * help of the MatrixAdapter and MultiVecAdapter classes) of the ordinal and
 * scalar types that will be used elsewhere in the concrete solver code.
 */

/**
 * \brief Amesos2::Solver: A templated interface for interaction with
 *        third-party direct sparse solvers.
 *
 * The Amesos2::Solver class is the statically polymorphic parent class of
 * each Amesos2 class named Amesos2_<i>SolverName</i> which interfaces with
 * the third-party solver named <i>SolverName</i>.
 *
 * <h1>Usage Examples</h1>
 *
 * ...
 */
template <template <class,class> class ConcreteSolver,
          class Matrix,
          class Vector >
class Solver : public  Amesos::SolverBase
{
public:

  // Grant public access to contained types
  typedef ConcreteSolver<Matrix,Vector>                            solver_type;
  typedef Matrix                                                   matrix_type;
  typedef Vector                                                   vector_type;
  typedef typename MatrixAdapter<matrix_type>::scalar_type         scalar_type;
  typedef typename MatrixAdapter<matrix_type>::local_ordinal_type  local_ordinal_type;
  typedef typename MatrixAdapter<matrix_type>::global_ordinal_type global_ordinal_type;
  typedef typename MatrixAdapter<matrix_type>::global_size_type    global_size_type;
  typedef typename MatrixAdapter<matrix_type>::node_type           node_type;


  /// \name Constructor/Destructor methods
  //@{

  /** \brief Initialize a Solver instance.
   *
   * A single constructor is supported, which accepts Teuchos::RCP objects
   * that point to a matrix-esque \x A, and vector-esque LHS object \X and RHS
   * object \x B.  This is the only constructor used by Amesos::Factory.
   *
   * \throw std::invalid_argument The shape of the matrix \c A is not
   * supported by the underlying solver.
   */
  Solver( Teuchos::RCP<Matrix> A, Teuchos::RCP<Vector> X, Teuchos::RCP<Vector> B );


  /// Destructor
  ~Solver( );


  /// Do not allow copying of this Solver object
  Solver(const solver_type& rhs);


  /// Do not allow copying of this Solver by assignment
  SolverBase& operator=(const solver_type* rhs);

  //@} End Constructor/Descructor block

  /// \name Mathematical functions
  //@{

  /**
   * \brief Pre-orders the matrix A for minimal fill-in
   *
   * Rearranges the rows and columns of the matrix A to minimize the amount of
   * fill-in of the non-zero entries of the matrix.  Pre-ordering may or may
   * not be supported by the underlying solver.  If not supported, a call to
   * this method simply does nothing.
   *
   * \return a reference to \c this .
   * \sa symbolicFactorization(), numericFactorization(), and solve()
   */
  SolverBase& preOrdering();


  /**
   * \brief Performs symbolic factorization on the matrix A.
   *
   * In addition to performing symbolic factorization on the matrix A, the call
   * to symbolicFactorization() implies that no change will be made to the
   * non-zero structure of the underlying matrix without a subsequent call to
   * symbolicFactorization().
   *
   * \pre
   *
   * \post
   * - Symbolic factorization will be performed (or marked to be performed)
   *   allowing numericFactorization() and solve() to be called.
   *
   * \return a reference to \c this .
   * \sa preOrdering(), numericFactorization(), and solve()
   */
  SolverBase& symbolicFactorization();


  /**
   * \brief Performs numeric factorization on the matrix A.
   *
   * In addition to performing numeric factorization on the matrix A, the call
   * to numericFactorization() implies that no change will be made to the
   * underlying matrix values without a subsequent call to
   * numericFactorization().
   *
   * \pre
   * - The non-zero structure of the matrix should not have changed since the
   *   last call to symbolicFactorization(). Other changes can have arbitrary
   *   consequences.
   * - The distribution of the matrix should not have changed
   *   since the last call to symbolicFactorization().
   *
   * \post
   *   Numeric factorization will be performed (or marked to be performed)
   *   allowing solve() to be performed correctly despite a potential change
   *   in the matrix values (though not in the non-zero structure).
   *
   * \return a reference to \c this
   * \sa preOrdering(), symbolicFactorization(), and solve()
   */
  SolverBase& numericFactorization();


  /**
   * \brief Solves \f$ A X = B\f$ (or \f$ A^T X = B\f$ )
   *
   * \pre
   * - The non-zero structure of the matrix should not have changed
   *   since the last call to symbolicFactorization().
   * - The distribution of the matrix should not have changed since
   *   the last call to symbolicFactorization().
   * - The matrix should not have changed since the last call to
   *   numericFactorization().
   *
   * \post
   *   X will be set such that \f$ A X = B\f$ (or \f$ A^T X = B\f$ ),
   *   within the limits of the accuracy of the underlying solver.
   *
   * \return void
   *
   * \sa preOrdering(), symbolicFactorization(), and numericFactorization()
   */
  void solve();


  /**
   * \brief Returns \c true if the solver can handle this matrix shape.
   *
   * Returns true if the matrix shape is one that the underlying
   * concrete sparse direct solver can handle. Classes that work
   * only on square matrices should return false for rectangular
   * matrices.  Classes that work only on symmetric matrices would
   * return false for non-symmetric matrices. etc.
   */
  bool matrixShapeOK();

  //@}  End Mathematical Functions group


  /// \name Parameter methods
  //@{

  /**
   * \brief Set/update internal variables and solver options.
   *
   * The setParameters method is consistent over all concrete solvers.  It
   * accepts general status and control parameters, as well as parameters
   * specific to a given solver.  If the solver does not recognize the
   * parameter, then it will simply be ignored
   *
   * \post
   * - Internal variables controlling the factorization and solve will be
   *   updated and take effect on all subsequent calls to
   *   numericFactorization() and solve().
   * - All parameters whose value is to differ from the default values must be
   *   included in \p parameterList. Parameters not specified in \p
   *   parameterList revert to their default values.
   *
   * \return a reference to \c this
   */
  SolverBase& setParameters(
    const Teuchos::RCP<Teuchos::ParameterList> & parameterList );

  /**
   * \example SimpleSolve_WithParameters.cpp
   *
   * An example of how to set solver parameters with the Superlu interface.
   */


  /**
   * \brief Return a const parameter list of all of the valid parameters that
   * this->setParameterList(...)  will accept.
   *
   * \note Check the documentation for your concrete solver to see a complete
   * list of the values that each parameter may take.  A solver may also
   * recognize multiple data types as arguments for a particular parameters
   * (eg. recognizing \c "YES" and \c "NO" as well as \c true and \c false ).
   */
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;


  /**
   * \brief Set or update internal variables and solver options.
   *
   * Redefined from Teuchos::ParameterListAcceptor
   *
   * \note Alias for \c setParameters()
   *
   * \param [in] parameterList
   */
  void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & parameterList)
    {
      setParameters(parameterList);
    }


  /**
   * \brief This is a empty stub
   *
   * \return
   */
  Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList()
    {
      return Teuchos::null;
    }


  /**
   * \brief This is an empty stub
   *
   * \return
   */
  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList()
    {
      return Teuchos::null;
    }


  //@} End Parameter methods


  /// \name Accessor methods
  //@{

  /// Returns a pointer to the Tpetra::Comm communicator with this operator.
  Teuchos::RCP<const Teuchos::Comm<int> > getComm() const
    {
      return matrixA_->getComm();
    }


  /// Returns the number of symbolic factorizations performed by this object.
  inline int getNumSymbolicFact() const { return( status_.numSymbolicFact_ ); }


  /// Returns the number of numeric factorizations performed by this object.
  inline int getNumNumericFact() const { return( status_.numNumericFact_ ); }


  /// Returns the number of solves performed by this object.
  inline int getNumSolve() const { return( status_.numSolve_ ); }

  //@} End Accessor methods


  /// Returns a short description of this Solver
  std::string description() const;


  /// Prints the status information about the current solver with some level
  /// of verbosity
  void describe(
    Teuchos::FancyOStream &out,
    const Teuchos::EVerbosityLevel verbLevel) const;


  /**
   * \brief Prints timing information about the current solver.
   *
   * The \c Amesos::Solver base class takes care of tracking total time spent
   * in the Amesos2 interface.  Concrete solver interface class are
   * responsible for reporting other timing statistics, which include time
   * spent in:
   *
   * - Redistribution of matrix objects,
   * - Conversion of matrix objects to solver-specific formats,
   * - Redistribution of multi-vector objects,
   * - Conversion of multi-vector objects to solver formats,
   * - TPL symbolic factorizations,
   * - TPL numeric factorizations, and
   * - TPL solves
   */
  void printTiming(
    Teuchos::FancyOStream &out,
    const Teuchos::EVerbosityLevel verbLevel) const;


  /**
   * \brief Extracts timing information from the current solver
   *
   * Results are placed into the parameter list \c timingParameterList.
   *
   * \param [out] timingParameterList Accepts timing information from the
   * current solver
   */
  void getTiming(Teuchos::ParameterList& timingParameterList) const;


  /**
   * \brief Return the name of this solver.
   *
   * The name is given by the underlying concrete solver instance
   *
   * \return A \c std::string which is the name of this solver
   */
  std::string name() const;


protected:

  /// The LHS operator
  Teuchos::RCP<MatrixAdapter<Matrix> >   matrixA_;

  /// The LHS vector/multi-vector
  Teuchos::RCP<MultiVecAdapter<Vector> > multiVecX_;

  /// The RHS vector/multi-vector
  Teuchos::RCP<MultiVecAdapter<Vector> > multiVecB_;

  /// Number of global rows in \c matrixA_
  global_size_type globalNumRows_;

  /// Number of global columns in \c matrixA_
  global_size_type globalNumCols_;

  /// Number of global non-zero values in \c matrixA_
  global_size_type globalNumNonZeros_;

  /* Status and Control data are handled in the Amesos::Status and
   * Amesos::Control base classes
   */

  /* Timers are handled in Amesos::Timers */

};				// End class Amesos::Solver


} // end namespace Amesos

#endif	// AMESOS2_SOLVER_DECL_HPP
