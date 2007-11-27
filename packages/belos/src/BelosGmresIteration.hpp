// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
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

#ifndef BELOS_GMRES_ITERATION_HPP
#define BELOS_GMRES_ITERATION_HPP

/*! \file BelosGmresIteration.hpp
    \brief Pure virtual base class which augments the basic interface for a Gmres linear solver iteration.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"
#include "BelosIteration.hpp"

namespace Belos {

  //! @name GmresIteration Structures 
  //@{ 
  
  /** \brief Structure to contain pointers to GmresIteration state variables.
   *
   * This struct is utilized by GmresIteration::initialize() and GmresIteration::getState().
   */
  template <class ScalarType, class MV>
  struct GmresIterationState {
    /*! \brief The current dimension of the reduction.
     *
     * This should always be equal to GmresIteration::getCurSubspaceDim()
     */
    int curDim;

    /*! \brief The current Krylov basis. */
    Teuchos::RCP<const MV> V;

    /*! \brief The current preconditioned Krylov basis (only used in flexible GMRES). */
    Teuchos::RCP<const MV> Z;

    /*! \brief The current Hessenberg matrix. 
     *
     * The \c curDim by \c curDim leading submatrix of H is the 
     * projection of problem->getOperator() by the first \c curDim vectors in V. 
     */

    Teuchos::RCP<const Teuchos::SerialDenseMatrix<int,ScalarType> > H;
    /*! \brief The current upper-triangular matrix from the QR reduction of H. */

    Teuchos::RCP<const Teuchos::SerialDenseMatrix<int,ScalarType> > R;
    /*! \brief The current right-hand side of the least squares system RY = Z. */

    Teuchos::RCP<const Teuchos::SerialDenseMatrix<int,ScalarType> > z;
    

    GmresIterationState() : curDim(0), V(Teuchos::null), Z(Teuchos::null),
			    H(Teuchos::null), R(Teuchos::null), 
			    z(Teuchos::null)
    {}
  };

  //@}

  //! @name GmresIteration Exceptions
  //@{ 
  
  /** \brief GmresIterationInitFailure is thrown when the GmresIteration object is unable to
   * generate an initial iterate in the GmresIteration::initialize() routine. 
   *
   * This std::exception is thrown from the GmresIteration::initialize() method, which is
   * called by the user or from the GmresIteration::iterate() method if isInitialized()
   * == \c false.
   *
   * In the case that this std::exception is thrown, 
   * GmresIteration::isInitialized() will be \c false and the user will need to provide
   * a new initial iterate to the iteration.
   */
  class GmresIterationInitFailure : public BelosError {public:
    GmresIterationInitFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};

  /** \brief GmresIterationOrthoFailure is thrown when the GmresIteration object is unable to
   * compute independent direction vectors in the GmresIteration::iterate() routine. 
   *
   * This std::exception is thrown from the GmresIteration::iterate() method.
   *
   */
  class GmresIterationOrthoFailure : public BelosError {public:
    GmresIterationOrthoFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};

  /** \brief GmresIterationLAPACKFailure is thrown when a nonzero return value is passed back
   * from an LAPACK routine.
   *
   * This std::exception is thrown from the GmresIteration::iterate() method.
   *
   */
  class GmresIterationLAPACKFailure : public BelosError {public:
    GmresIterationLAPACKFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};
  
  //@}


template<class ScalarType, class MV, class OP>
class GmresIteration : virtual public Iteration<ScalarType,MV,OP> {

  public:

  //! @name State methods
  //@{ 
  /*! \brief Initialize the solver to an iterate, providing a complete state.
   *
   * The %GmresIteration contains a certain amount of state, consisting of the current 
   * residual, preconditioned residual, and decent direction.
   *
   * initialize() gives the user the opportunity to manually set these,
   * although only the current unpreconditioned residual is required.
   *
   * \post 
   * <li>isInitialized() == \c true (see post-conditions of isInitialize())
   *
   * \note For any pointer in \c newstate which directly points to the multivectors in 
   * the solver, the data is not copied.
   */
  virtual void initializeGmres(GmresIterationState<ScalarType,MV> newstate) = 0;

  /*! \brief Get the current state of the linear solver.
   *
   * The data is only valid if isInitialized() == \c true.
   *
   * \returns A GmresIterationState object containing const pointers to the current solver state.
   */
  virtual GmresIterationState<ScalarType,MV> getState() const = 0;
  //@}

  //! @name Status methods
  //@{ 

  //! Method for updating QR factorization of upper Hessenberg matrix
  /*! \note If \c dim >= \c getCurSubspaceDim() and \c dim < \c getMaxSubspaceDim(), then 
            the \c dim-th equations of the least squares problem will be updated.
  */
  virtual void updateLSQR( int dim = -1 ) = 0;

  //! Get the dimension of the search subspace used to generate the current solution to the linear problem.
  virtual int getCurSubspaceDim() const = 0; 

  //! Get the maximum dimension allocated for the search subspace.
  virtual int getMaxSubspaceDim() const = 0;

  //@}

  //! @name Accessor methods
  //@{ 

  /*! \brief Set the blocksize and number of blocks to be used by the
   * iterative solver in solving this linear problem.
   *
   *  Changing either the block size or the number of blocks will reset the
   *  solver to an uninitialized state.
   */
  virtual void setSize(int blockSize, int numBlocks) = 0;

  //@}

};

} // end Belos namespace

#endif /* BELOS_GMRES_ITERATION_HPP */
