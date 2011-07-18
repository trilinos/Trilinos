// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2010 Sandia Corporation
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

/**
  \file   Amesos2_Solver_decl.hpp
  \author Eric T Bavier <etbavier@sandia.gov>
  \date   Thu Dec 31 10:32:08 2009

  \brief  Interface to Amesos2 solver objects
*/

#ifndef AMESOS2_SOLVER_DECL_HPP
#define AMESOS2_SOLVER_DECL_HPP

#include <Teuchos_Describable.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>

#include "Amesos2_Status.hpp"

namespace Amesos2 {


  /**
   * \brief Interface to Amesos2 solver objects
   *
   * Specifies a uniform interface for interaction with Amesos2 solver wrappers
   * to third-party libraries.
   *
   * This class holds no definitions itself; it is a contract for
   * concrete solver interfaces only.
   *
   * \ingroup amesos2_solver_framework
   */
  template <class Matrix, class Vector>
  class Solver : public Teuchos::Describable {

  public:

    typedef Solver<Matrix,Vector> type;

    /// \name Mathematical Functions
    //@{

    /** \brief Pre-orders the matrix.
     *
     * Uses the default solver option, unless a solver-specific
     * pre-ordering parameter is given.
     *
     * \sa setParameters
     */
    virtual type& preOrdering( void ) = 0;


    /** \brief Performs symbolic factorization on the matrix
     *
     * \pre
     *  - The matrix A must not be \c null
     */
    virtual type& symbolicFactorization( void ) = 0;


    /** \brief Performs numeric factorization on the matrix
     *
     * numericFactorization checks first that symbolicFactorization
     * has successfully been called, and if not, calls it before
     * continuing.
     * 
     * \pre
     *  - The matrix A must not be \c null
     *
     * \post
     *  - The factors L and U of A are computed
     */
    virtual type& numericFactorization( void ) = 0;


    /** \brief Solves \f$ A X = B\f$ (or \f$ A^T X = B\f$ )
     *
     * solve checks first that numericFactorization has successfully
     * been called, and if not, calls it before continuing.
     * 
     * \pre
     *  - The (multi)vectors \c X and \c B must not be \c null
     *
     * \post
     *  - The (multi)vector \c X (given at construction time) contains
     *    the solution to the system.
     */
    virtual void solve( void ) = 0;


    /** \brief Solve \f$ A X = B\f$ using the given X and B vectors.
     *
     * This overload of solve uses the given X and B vectors when
     * solving.  This X and B are used in place of any X and B that
     * were given upon construction of the Amesos2 solver instance and
     * are used only for this solve.
     *
     * If a permanent change of X and B are required, see the setX()
     * and setB() methods.
     *
     * \post
     *  - The (multi)vector \c X contains the solution to the system
     *  - The \c X and \c B given at construction time (if any) are unchanged.
     */
    virtual void solve(const Teuchos::Ptr<Vector>       X,
		       const Teuchos::Ptr<const Vector> B) const = 0;


    /** \brief Solve \f$ A X = B\f$ using the given X and B vectors.
     *
     * This overload of solve uses the given X and B vectors when
     * solving.  This X and B are used in place of any X and B that
     * were given upon construction of the Amesos2 solver instance and
     * are used only for this solve.
     *
     * If a permanent change of X and B are required, see the setX()
     * and setB() methods.
     *
     * \post
     *  - The (multi)vector \c X contains the solution to the system
     *  - The \c X and \c B given at construction time (if any) are unchanged.
     */
    virtual void solve(Vector* X, const Vector* B) const = 0;

    //@} End Mathematical Functions


    /** \name Parameter Methods
     * @{
     */

    /** \brief Set/update internal variables and solver options.
     *
     * Expects that parameterList be named "Amesos2".  That list may
     * contain Amesos2-specific parameters.  In addition, it may
     * contain sublist for solver-specific parameters.  These sublists
     * should be named according to what is returned by the name()
     * function (i.e. The solver's name when enabling for Amesos2
     * during configuration).
     *
     * See each solver interface directly for a list of the supported
     * parameters for that solver.
     */
    virtual type& setParameters( const Teuchos::RCP<Teuchos::ParameterList> & parameterList ) = 0;


    /**
     * \brief Return a const parameter list of all of the valid parameters that
     * this->setParameterList(...)  will accept.
     */
    virtual Teuchos::RCP<const Teuchos::ParameterList> getValidParameters( void ) const = 0;

    /// @} End Parameter Methods


    /** \name Accessor Methods
     * @{
     */

    /** \brief Sets the matrix A of this solver
     *
     * Calling this method will reset any internal state of the solver
     * which was dependent on the previous matrix.
     *
     * Passing Teuchos::null is a special case, however.  This will
     * release the matrix, but will not reset any internal state.  It
     * is expected then, that symbolicFactorization or
     * numericFactorization will not be called again until another
     * matrix is set using this method.  This allows the user to free
     * a matrix from memory once the L and U factors have been
     * computed.  For example:
     *
     * \code
     * Solver<MAT,MV> solver = create(A);
     * solver->symbolicFactorization().numericFactorization();
     * A = Teuchos::null;
     * solver->setA(Teuchos::null); // release A
     * solver->solve(X, B);
     * \endcode
     */
    virtual void setA( const Teuchos::RCP<const Matrix> a ) = 0;


    virtual void setA( const Matrix* a ) = 0;
    

    /// Returns \c true if the solver can handle the matrix shape
    virtual bool matrixShapeOK( void ) = 0;


    /// Sets the LHS vector X
    virtual void setX( const Teuchos::RCP<Vector> x ) = 0;


    /// Sets the LHS vector X using a raw pointer
    virtual void setX( Vector* x ) = 0;


    /// Returns the vector that is the LHS of the linear system
    virtual const Teuchos::RCP<Vector> getX( void ) = 0;


    /// Returns a raw pointer to the LHS of the linear system
    virtual Vector* getXRaw( void ) = 0;


    /// Sets the RHS vector B
    virtual void setB( const Teuchos::RCP<const Vector> b ) = 0;


    /// Sets the RHS vector B using a raw pointer
    virtual void setB( const Vector* b ) = 0;


    /// Returns the vector that is the RHS of the linear system
    virtual const Teuchos::RCP<const Vector> getB( void ) = 0;


    /// Returns a raw pointer to the RHS of the linear system
    virtual const Vector* getBRaw( void ) = 0;


    /// Returns a pointer to the Teuchos::Comm communicator with this matrix
    virtual Teuchos::RCP<const Teuchos::Comm<int> > getComm( void ) const = 0;


    /// Returns a reference to this solver's internal status object
    virtual Status& getStatus() const = 0;


    /// Returns the number of pre-orderings performed by this object.
    virtual int getNumPreOrder( void ) const = 0;
    

    /// Returns the number of symbolic factorizations performed by this object.
    virtual int getNumSymbolicFact( void ) const = 0;


    /// Returns the number of numeric factorizations performed by this object.
    virtual int getNumNumericFact( void ) const = 0;


    /// Returns the number of solves performed by this object.
    virtual int getNumSolve( void ) const = 0;


    /// Return the name of this solver.
    virtual std::string name( void ) const = 0;

    /// @} End Accessor Methods


    /** \name Methods implementing Describable
     * @{
     */

    /// Returns a short description of this Solver
    virtual std::string description( void ) const = 0;


    /// Prints the status information about the current solver with some level
    /// of verbosity.
    virtual void describe( Teuchos::FancyOStream &out,
			   const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default ) const = 0;

    /// @} End Methods implementing Describable


    /** \name Performance and Timing
     * @{
     */
    
    /// Prints timing information about the current solver.
    virtual void printTiming( Teuchos::FancyOStream &out,
			      const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default ) const = 0;


    /**
     * \brief Extracts timing information from the current solver.
     *
     * Results are placed in the parameter list \c timingParameterList
     *
     * \param timingParameterList Accepts timing information from the
     * current solver
     */
    virtual void getTiming( Teuchos::ParameterList& timingParameterList ) const = 0;

    /// @} End Performance and Timing

  };                              // End class Solver


} // end namespace Amesos2

#endif	// AMESOS2_SOLVER_BASE_DECL_HPP
