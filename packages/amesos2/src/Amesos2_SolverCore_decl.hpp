// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package 
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

/**
  \file   Amesos2_SolverCore_decl.hpp
  \author Eric T Bavier <etbavier@sandia.gov>
  \date   Thu May 27 13:51:05 CDT 2010

  \brief  Templated core-functionality class for Amesos2 solvers
*/

#ifndef AMESOS2_SOLVERCORE_DECL_HPP
#define AMESOS2_SOLVERCORE_DECL_HPP

#include <string>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "Amesos2_Solver.hpp"
#include "Amesos2_MatrixTraits.hpp"
#include "Amesos2_MatrixAdapter_decl.hpp"
#include "Amesos2_MultiVecAdapter_decl.hpp"
#include "Amesos2_TypeDecl.hpp"

#include "Amesos2_Control.hpp"
#include "Amesos2_Status.hpp"
#include "Amesos2_Timers.hpp"

namespace Amesos2 {


  /* This is the base class to be used in a *statically* polymorphic
   * way. E.g. for the Superlu solver:
   *
   * In Amesos2_Superlu.hpp:
   * class Superlu : SolverCore<Superlu> { ... }
   *
   * Each concrete solver will implement several private sub-functions
   * that will be called within the common code for each function.
   */

  /**
   * \brief Amesos2::SolverCore: A templated interface for interaction
   *        with third-party direct sparse solvers.
   *
   * The Amesos2::SolverCore class is the statically polymorphic parent
   * class of each Amesos2 class named Amesos2_<i>SolverName</i> which
   * interfaces with the third-party solver named <i>SolverName</i>.
   *
   * Each concrete solver interface implements several private
   * "sub-functions" that are called within common code for that
   * function.
   *
   * This static base class also provides a convenient way to handle the
   * template parameters and private typedef'ing (with the help of the
   * Matrix and MultiVector adapter classes) of the ordinal and scalar
   * types that can be used elsewhere in the concrete solver code.
   *
   * \ingroup amesos2_solver_framework
   */
  template <template <class,class> class ConcreteSolver,
	    class Matrix,
	    class Vector >
  class SolverCore : public Amesos2::Solver<Matrix,Vector>
  {
  public:

    // Grant public access to contained types
    typedef SolverCore<ConcreteSolver,Matrix,Vector>                       type;
    typedef Solver<Matrix,Vector>                                    super_type;
    typedef ConcreteSolver<Matrix,Vector>                           solver_type;
    typedef Matrix                                                  matrix_type;
    typedef Vector                                                  vector_type;
    typedef typename MatrixAdapter<matrix_type>::scalar_t           scalar_type;
    typedef typename MatrixAdapter<matrix_type>::local_ordinal_t    local_ordinal_type;
    typedef typename MatrixAdapter<matrix_type>::global_ordinal_t   global_ordinal_type;
    typedef typename MatrixAdapter<matrix_type>::global_size_t      global_size_type;
    typedef typename MatrixAdapter<matrix_type>::node_t             node_type;


    /// \name Constructor/Destructor methods
    //@{

    /** \brief Initialize a Solver instance.
     *
     * A single constructor is supported, which accepts Teuchos::RCP
     * objects that point to a matrix-esque \x A, and vector-esque LHS
     * object \X and RHS object \x B.  This is the only constructor used
     * by Amesos2::create.
     *
     * \throw std::invalid_argument The shape of the matrix \c A is not
     * supported by the underlying solver.
     */
    SolverCore( Teuchos::RCP<const Matrix> A,
		Teuchos::RCP<Vector>       X,
		Teuchos::RCP<const Vector> B );


    /// Destructor
    ~SolverCore( );


    /// Do not allow copying of this Solver object
    SolverCore(const solver_type& rhs);


    /// Do not allow copying of this Solver by assignment
    super_type& operator=(const solver_type* rhs);

    //@} End Constructor/Destructor block

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
    super_type& preOrdering();


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
    super_type& symbolicFactorization();


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
    super_type& numericFactorization();


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


    void solve(const Teuchos::Ptr<Vector> X, const Teuchos::Ptr<const Vector> B) const;


    void solve(Vector* X, const Vector* B) const;
  

    //@}  End Mathematical Functions group


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

    void setA( const Teuchos::RCP<const Matrix> a, EPhase keep_phase = CLEAN );

    void setA( const Matrix* a, EPhase keep_phase = CLEAN ){ setA(Teuchos::rcp(a), keep_phase); }

    void setX(const Teuchos::RCP<Vector> x){ multiVecX_ = x; }

    void setX(Vector* x){ multiVecX_ = Teuchos::rcp(x); }
  
    const Teuchos::RCP<Vector> getX(){ return( multiVecX_ ); }

    Vector* getXRaw(){ return multiVecX_.getRawPtr(); }

    void setB(const Teuchos::RCP<const Vector> b){ multiVecB_ = b; }

    void setB(const Vector* b){ multiVecB_ = Teuchos::rcp(b); }

    const Teuchos::RCP<const Vector> getB(){ return( multiVecB_ ); }

    const Vector* getBRaw(){ return multiVecB_.getRawPtr(); }


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
    super_type& setParameters(
			      const Teuchos::RCP<Teuchos::ParameterList> & parameterList );


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

    /// Returns a pointer to the Teuchos::Comm communicator with this operator.
    Teuchos::RCP<const Teuchos::Comm<int> > getComm() const
    {
      return matrixA_->getComm();
    }

  
    /// Returns a reference to this solver's internal status object
    inline Status& getStatus() const { return( status_ ); }


    //@} End Accessor methods


    /// Returns a short description of this Solver
    std::string description() const;


    /// Prints the status information about the current solver with some level
    /// of verbosity
    void describe(Teuchos::FancyOStream &out,
		  const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;


    /**
     * \brief Prints timing information about the current solver.
     *
     * The \c Amesos2::SolverCore base class takes care of tracking
     * total time spent in the Amesos2 interface.  Concrete solver
     * interface class are responsible for reporting other timing
     * statistics, which include time spent in:
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

  private:

    /** \brief Refresh this solver's internal data about A
     * 
     * Called whenever it would be necessary to refresh a solver's
     * internal storage of the matrix A, which is whenever a phase is
     * called that is equal to or below the current call.
     *
     * For example, say a user has just previously called solve(), then
     * calls numericFactorization().  Since the solve phase is greater
     * than the numericFactorization phase, this is an indication that
     * the internal store of A needs refreshing, since the user
     * (assuming the user know what she's doing) wouldn't otherwise need
     * to call numericFactorization following a solve.
     */
    void loadA(EPhase current_phase);

  protected:

    /** \brief Set the number of non-zero values in the \f$L\f$ and \f$U\f$ factors
     *
     * Concrete solver classes may call this method if they wish to
     * (or are able to) report the number of conbined non-zero count
     * for the \f$L\f$ and \f$U\f$ factors.
     */
    void setNnzLU(size_t nnz){ status_.lu_nnz_ = nnz; }

    /// The LHS operator
    Teuchos::RCP<const MatrixAdapter<Matrix> > matrixA_;

    /** \internal
     * 
     * If \c true indicates that the current matrix A has been loaded
     * into internal solver structures.
     */
    bool matrix_loaded_;


    /// The LHS vector/multi-vector
    Teuchos::RCP<Vector> multiVecX_;

    /** \brief The RHS vector/multi-vector
     * 
     * We point to a const Vector because Amesos2 should never directly
     * modify B.
     */
    Teuchos::RCP<const Vector> multiVecB_;

    /// Number of global rows in \c matrixA_
    global_size_type globalNumRows_;

    /// Number of global columns in \c matrixA_
    global_size_type globalNumCols_;

    /// Number of global non-zero values in \c matrixA_
    global_size_type globalNumNonZeros_;

    /// Index base of rowmap of \c matrixA_
    global_size_type rowIndexBase_;

    /// Index base of column map of \c matrixA_
    global_size_type columnIndexBase_;


    /// Holds status information about a solver
    mutable Status status_;

    /// Parameters for solving
    Control control_;

    /// Various timing statistics
    mutable Timers timers_;


    /* Useful MPI vars */

    /// The MPI rank of this image
    int rank_;

    /// If \c true, then this is the root processor
    bool root_;

    /// Number of process images in the matrix communicator
    int nprocs_;

  };				// End class Amesos2::SolverCore


} // end namespace Amesos2

#endif	// AMESOS2_SOLVERCORE_DECL_HPP
