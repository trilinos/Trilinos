// @HEADER
// ***********************************************************************
//
//                Amesos: Direct Sparse Solver Package
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
#include <Teuchos_ParameterListAcceptor.hpp>

#include <Tpetra_Operator.hpp>
#include <Tpetra_RowMatrix.hpp>
#include <Tpetra_MultiVector.hpp>

#include "Amesos2_SolverBase.hpp"
#include "Amesos2_MatrixAdapter.hpp"
#include "Amesos2_MultiVecAdapter.hpp"

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
/** \brief Amesos2::Solver: A templated interface for interaction with
 *         third-party direct sparse solvers.
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
    
  /* Constructors will accept either pointers to matrices and vectors or
   * Teuchos::RCP to matrices and vectors, although only the latter is
   * returned by Amesos::Factory
   */
  Solver( void );	

  Solver( Teuchos::RCP<Matrix> A, Teuchos::RCP<Vector> X, Teuchos::RCP<Vector> B );

  /// Destructor
  ~Solver( );

  /// Do not allow copying of this Solver object
  Solver(const solver_type& rhs);

  /// Do not allow copying of this Solver by assignment
  SolverBase& operator=(const solver_type* rhs);

  //@{ \name Mathematical functions

  /** \brief Pre-orders the matrix A for minimal fill-in
   * 
   * Rearranges the rows and columns of the matrix A to minimize the amount of
   * fill-in of the non-zero entries of the matrix.  Pre-ordering may or may
   * not be supported by the underlying solver.  If not supported, a call to
   * this method simply does nothing.
   *
   * \return a reference to \c this.
   * \sa symbolicFactorization(), numericFactorization(), and solve()
   */
  SolverBase& preOrdering();
    

  /** \brief Performs symbolic factorization on the matrix A.
   *
   * In addition to performing symbolic factorization on the matrix A, the call
   * to symbolicFactorization() implies that no change will be made to the
   * non-zero structure of the underlying matrix without a subsequent call to
   * symbolicFactorization().
   *  
   * \pre
   * - matrixA != null (return -1)
   * - matrixShapeOK() == true (return -6)
   * 
   * \post
   * - Symbolic factorization will be performed (or marked to be
   * performed) allowing numericFactorization() and solve() to
   * be called.
   *     
   * \return a reference to \c this
   * \sa preOrdering(), numericFactorization(), and solve()
   */
  SolverBase& symbolicFactorization();
  

  /** \brief Performs numeric factorization on the matrix A.
   * In addition to performing numeric factorization on the matrix
   * A, the call to numericFactorization() implies that no change
   * will be made to the underlying matrix without a subsequent call
   * to numericFactorization().
   *
   * \pre
   * - matrixA != 0 (return -1)
   * - matrixShapeOK() == true (return -6)
   * - The non-zero structure of the matrix should not have changed since the
   *   last call to symbolicFactorization(). (return -2 if the number of
   *   non-zeros changes).  Other changed can have arbitrary consequences.
   * - The distribution of the matrix should not have changed
   *   since the last call to symbolicFactorization().
   * - The matrix should be indexed from 0 to n-1, unless the
   *   parameter "Reindex" was set to "true" prior to the call to
   *   symbolicFactorization(). (return -3 if caught)
   * - The parameter "Reindex" should not be set to "true" except
   *   on CrsMatrices. (return -4)
   * ->check	 - The parameter "Reindex" should not be set to "true" unless
   *             Amesos was built with EpetraExt, i.e. with
   *             --enable-epetraext on the configure line. (return -4)
   * - Internal errors return -5
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


  /** \brief Solves \f$ A X = B\f$ (or \f$ A^T X = B\f$ ) 
   * \pre
   * - matrixA != null (return -1)
   * - matrixShapeOK() == true (return -6)
   * - The non-zero structure of the matrix should not have changed
   *   since the last call to symbolicFactorization().
   * - The distribution of the matrix should not have changed since
   *   the last call to symbolicFactorization.
   * - The matrix should not have changed since the last call to
   *   numericFactorization().
   * 
   * \post
   *   X will be set such that \f$ A X = B\f$ (or \f$ A^T X = B\f$ ),
   *   within the limits of the accuracy of the underlying solver.
   *   
   * \return void
   * \sa preOrdering(), symbolicFactorization(), and numericFactorization()
   */
  void solve();

  
  /** \brief Returns true if the solver can handle this matrix shape.
   *
   * Returns true if the matrix shape is one that the underlying
   * concrete sparse direct solver can handle. Classes that work
   * only on square matrices should return false for rectangular
   * matrices.  Classes that work only on symmetric matrices would
   * return false for non-symmetric matrices.
   */
  bool matrixShapeOK();
  

  //@}  End Mathematical Functions group


  //@{ \name Attribute set methods

  /** \brief Updates internal variables.
   *
   * \pre
   * None
   * 
   * \post
   * - Internal variables controlling the factorization and solve
   * will be updated and take effect on all subsequent calls to
   * numericFactorization() and solve().
   * - All parameters whose value is to differ from the default
   * values must be included in \p parameterList. Parameters not
   * specified in \p parameterList revert to their default values.
   * 
   * \param parameterList 
   *   
   * \return a reference to \c this
   */
  SolverBase& setParameters( const Teuchos::RCP<Teuchos::ParameterList> & parameterList );
  

  /// Returns a pointer to the Tpetra::Comm communicator with this operator.
  inline const Teuchos::RCP<const Teuchos::Comm<int> >& getComm() const;
  

  /// Returns the number of symbolic factorizations performed by this object.
  inline int getNumSymbolicFact() const { return(status_.numSymbolicFact_); }


  /// Returns the number of numeric factorizations performed by this object.
  inline int getNumNumericFact() const { return(status_.numNumericFact_); }


  /// Returns the number of solves performed by this object.
  inline int getNumSolve() const { return(status_.numSolve_); }


  /// Returns a short description of this Solver
  std::string description() const;

  
  /// Prints the status information about the current solver with some level
  /// of verbosity
  void describe(
    Teuchos::FancyOStream &out,
    const Teuchos::EVerbosityLevel verbLevel) const;


  /// Prints timing information about the current solver.
  void printTiming(
    Teuchos::FancyOStream &out,
    const Teuchos::EVerbosityLevel verbLevel) const;

  
  /** 
   * Redefined from Teuchos::ParameterListAcceptor
   *  
   * \param [in] parameterList 
   */
  void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & parameterList)
    {
      setParameters(parameterList);
    }
  

  /** 
   * This is a empty stub
   *   
   * \return 
   */
  Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList()
    {
      return Teuchos::null;
    }
  

  /** 
   * This is an empty stub
   *    
   * \return 
   */
  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList()
    {
      return Teuchos::null;
    }
  

  /** \brief Extracts timing information from the current solver
   *
   * Results are placed into the parameter list \c timingParameterList.
   *   
   * \param [out] timingParameterList Accepts timing information from the
   * current solver 
   */
  void getTiming(Teuchos::ParameterList& timingParameterList) const;

 
  /** \brief Return the name of this solver.
   * 
   * The name is given by the underlying concrete solver instance
   * 
   * \return A \c std::string which is the name of this solver
   */
  std::string name() const;
  

protected:

  Teuchos::RCP<MatrixAdapter<Matrix> >   matrixA_;
  Teuchos::RCP<MultiVecAdapter<Vector> > multiVecX_;
  Teuchos::RCP<MultiVecAdapter<Vector> > multiVecB_;

  global_size_type globalNumRows_;
  global_size_type globalNumCols_;
  global_size_type globalNumNonZeros_;

  // Status and Control data are handled in the Amesos2::Status and
  // Amesos2::Control base classes

  // Timers are handled in Amesos2::Timers
};				// End class Amesos::Solver

  
} // end namespace Amesos

#endif	// AMESOS2_SOLVER_DECL_HPP
