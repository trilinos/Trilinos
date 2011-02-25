//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
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
// ************************************************************************
//@HEADER

#ifndef __Belos_GmresBase_hpp
#define __Belos_GmresBase_hpp

/// \file BelosGmresBase.hpp 
/// \brief Common state and functionality for new GMRES implementations
/// \author Mark Hoemmen

#include <BelosLinearProblem.hpp>
#include <BelosOrthoManager.hpp>
#include <BelosMultiVecTraits.hpp>

#include <Teuchos_BLAS.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_SerialDenseVector.hpp>

#include <algorithm>
#include <iterator>
#include <utility> // std::pair

namespace Belos {

  /// \class GmresCantExtendBasis
  /// \author Mark Hoemmen
  /// \brief Raised by GmresBase::extendBasis() if can't extend basis
  ///
  /// Implementations of GmresBase's extendBasis() method must raise
  /// this exception if no more candidate basis vectors can be
  /// generated.  The cause should not be due to the numerical values
  /// in the candidate basis vector(s); the orthogonalization will
  /// detect that.  The usual cause is that the allotted maximum
  /// number of basis vectors has been reached.  Subclasses may choose
  /// instead to increase this maximum number and attempt to
  /// reallocate storage.
  ///
  /// \note BelosError is a subclass of std::logic_error.
  /// GmresCantExtendBasis "is a" logic_error, because callers of
  /// GmresBase::advance() should use the canAdvance() method rather
  /// than a try/catch to limit the number of iterations.  GmresBase
  /// is an implementation class (with the interface that I want).  It
  /// will be wrapped by a subclass of Belos::Iteration (with the
  /// interface I don't want but have to use for the sake of backwards
  /// compatibility).  The Iteration subclass will use the status test
  /// to control iteration and limit the number of iterations to the
  /// maximum number accepted by GmresBase.  So, GmresCantExtendBasis
  /// should never be thrown by correct code, thus it's a logic_error.
  ///
  /// \note Incidentally, Belos' exception hierarchy is a bit broken,
  /// since BelosError inherits from std::logic_error.  Logic errors
  /// are programmer bugs, and thus attempting to recover from them is
  /// a bad idea.  However, one can imagine some iterative method
  /// errors from which recovery is possible.  These should not be
  /// considered logic errors.
  ///
  class GmresCantExtendBasis : public BelosError {
  public:
    GmresCantExtendBasis(const std::string& what_arg) : BelosError(what_arg) {}
  };

  /// \brief Candidate "basis" isn't a basis
  ///
  /// Thrown by GmresBase::advance(), when it rejects the computed
  /// candidate basis vector(s) due to (numerical) rank deficiency, 
  /// and doesn't know how to recover.
  ///
  /// This usually means that after orthogonalizing the candidate
  /// basis vector(s) from extendBasis(), they are not full rank.  In
  /// the case of standard GMRES, this means the candidate basis
  /// vector has zero norm.  For CA-GMRES, the vectors might have
  /// nonzero norm, but are not full rank.  CA-GMRES may choose to
  /// retry with a shorter candidate basis length, but if the
  /// candidate basis length is too short, it may opt to "give up."
  /// In that case, advance() throws this exception.  Restarting with
  /// standard GMRES may be a good idea in that case.
  ///
  /// Applications may choose to recover from or deal with this error
  /// in one or more of the following ways: 
  ///
  /// - For CA-GMRES: equilibrate the matrix (and preconditioner)
  ///   and try again, or restart using standard GMRES
  /// - Assume that the (preconditioned) matrix is singular, 
  ///   and use an appropriate solver
  /// - Solve again with higher floating-point precision
  ///
  /// It might be good to verify that the matrix (and preconditioner)
  /// are nonzero.
  ///
  /// \note This may not necessarily be a "logic error" (i.e., coding
  ///   bug), since it can result from a matrix which is nonsingular
  ///   in exact arithmetic, but ill-conditioned.  However, this
  ///   exception must inherit from BelosError, which is an
  ///   std::logic_error.
  class GmresRejectsCandidateBasis : public BelosError {
  public:
    GmresRejectsCandidateBasis(const std::string& what_arg) : BelosError(what_arg) {}
  };


  /// \class GmresBase
  /// \brief Common state and functionality for new GMRES implementations
  /// \ingroup belos_solver_framework
  /// \author Mark Hoemmen
  ///
  /// This class includes both state and functionality that are useful
  /// for different implementations of the Generalized Minimal
  /// Residual (GMRES) method of Saad and Schultz, for iterative
  /// solution of nonsymmetric nonsingular linear systems.  Both the
  /// original GMRES algorithm (Saad and Schultz 1986) and Flexible
  /// GMRES (FGMRES) (Saad 1993) are supported.  This class does not
  /// implement the actual iterations; this is left for subclasses.
  /// Furthermore, it does not implement features like recycling, but
  /// it does include hooks for subclasses to add in such features.
  ///
  /// References:
  /// - Saad and Schultz, "GMRES: A generalized minimal residual
  ///   algorithm for solving nonsymmetric linear systems", SISSC,
  ///   vol. 7, pp. 856-869, 1986.
  /// - Saad, "A flexible inner-outer preconditioned GMRES algorithm",
  ///   SISC, vol. 14, pp. 461-469, 1993.
  ///
  template<class Scalar, class MV, class OP>
  class GmresBase {
  public:
    typedef Scalar scalar_type;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;
    typedef MV multivector_type;
    typedef OP operator_type;

  private:
    typedef Teuchos::SerialDenseMatrix<int,Scalar> mat_type;
    typedef Teuchos::SerialDenseVector<int,Scalar> vec_type;
    typedef Teuchos::BLAS<int, scalar_type> blas_type;
    typedef MultiVecTraits<scalar_type, MV> MVT;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef Teuchos::ScalarTraits<magnitude_type> STM;

  public:
    /// \brief Constructor
    ///
    /// Construct and initialize a new GmresBase iteration to solve a
    /// given linear problem.
    ///
    /// \param problem [in/out] Linear problem.  On input, we use the
    ///   starting guess (x0) and the initial residual (r0) to
    ///   initialize the iteration.  The iteration may call
    ///   updateSolution().  On output, if the solution has been
    ///   updated, the vector returned by getLHS() will be modified.
    /// \param ortho [in] Orthogonalization manager
    /// \param maxIterCount [in] Maximum number of iterations before
    ///   restart.  The number of vectors' worth of storage this
    ///   constructor allocates is proportional to this, so choose
    ///   carefully.
    /// \param flexible [in] Whether or not to run the Flexible
    ///   variant of GMRES (FGMRES).  This requires twice as much
    ///   vector storage, but lets the preconditioner change in every
    ///   iteration.  This only works with right preconditioning, not
    ///   left or split preconditioning.
    ///
    /// \note Instantiating a new GmresBase subclass instance is the
    /// only way to change the linear problem.  This ensures that the
    /// V (and Z, for Flexible GMRES) spaces (especially their data
    /// distributions) are correct.  It's allowed to change a
    /// LinearProblem instance's operator and right-hand, for example,
    /// and afterwards they might not even have the same number of
    /// rows, let alone the same data distribution as before.
    /// Belos::MultiVecTraits offers a limited interface to
    /// multivectors and operators that doesn't let us access the data
    /// distribution, so we choose instead to reinitialize from
    /// scratch whenever the linear problem is changed.  The way to
    /// avoid this reinitialization would be to refactor Belos to use
    /// only Thyra's MultiVectorBase and LinearOpBase, instead of
    /// being templated on MV and OP.  Thyra's interfaces reify and
    /// offer (abstract) access to the vector space (VectorSpaceBase)
    /// to which (multi)vectors belong, and the domain and range of an
    /// operator.
    GmresBase (const Teuchos::RCP<LinearProblem<Scalar, MV, OP> >& problem,
	       const Teuchos::RCP<const OrthoManager<Scalar, MV> >& ortho,
	       const int maxIterCount,
	       const bool flexible);

    //! Whether it is legal to call advance()
    bool canAdvance () const { return canExtendBasis(); }

    /// \brief Perform the next group of iteration(s)
    ///
    /// Perform the next group of iteration(s), without any
    /// convergence or termination tests.  Throw GmresCantExtendBasis
    /// if there is no more room for basis vectors; throw
    /// GmresRejectsCandidateBasis if the candidate basis vector(s) can't be
    /// accepted.
    ///
    /// \note Different subclasses may compute different numbers of
    /// iteration(s) in this method.  Standard (F)GMRES attempts
    /// exactly one iteration, whereas CA-(F)GMRES does work
    /// equivalent to several iterations of standard (F)GMRES.
    ///
    void advance();

    /// Accept the candidate basis and prepare for the next iteration.
    /// 
    /// \note Subclasses that produce more than one basis vector at a
    ///   time should override this implementation appropriately.
    virtual void acceptCandidateBasis () {
      curNumIters_++;
    }

    /// Reject the candidate basis
    /// 
    /// \note Subclasses that retry rejected iterations should
    ///   override this implementation appropriately.
    virtual void rejectCandidateBasis () {
      std::ostringstream os;
      os << (flexible_ ? "Flexible GMRES" : "GMRES") 
	 << " rejects the computed candidate basis vector(s)";
      TEST_FOR_EXCEPTION(true, GmresRejectsCandidateBasis, os.str());
    }

    /// \brief The number of completed iteration(s) (>= 0)
    /// 
    /// Does not include the starting vector.
    int getNumIters() const { return curNumIters_; }

    /// \brief Maximum number of iterations allowed before restart
    ///
    /// Does not include the starting vector.
    int maxNumIters() const { return maxNumIters_; }

    /// \brief Restart, and optionally reset max iteration count
    /// 
    /// First run preRestartHook().  Then, update the LinearProblem
    /// with the current approximate solution, and let it recompute
    /// the (exact, not "native") residual.  Restart, reallocating
    /// space for basis vectors and the projected least-squares
    /// problem if necessary.  Finally, run postRestartHook().  The
    /// hooks let subclasses implement things like subspace recycling.
    ///
    /// \param maxIterCount [in] New maximum iteration count for the
    ///   next restart cycle.  If same as the old iteration count (the
    ///   default), we reuse existing space for the Krylov basis
    ///   vectors.
    ///
    /// \note The restart will be invalid if the LinearProblem changes
    ///   between restart cycles (except perhaps for the
    ///   preconditioner if we are running Flexible GMRES).  It's
    ///   impossible for iteration classes to enforce this
    ///   requirement, since LinearProblem objects are mutable.  If
    ///   you want to solve a different linear problem (e.g., with a
    ///   different matrix A or right-hand side b), instantiate a new
    ///   subclass of GmresBase.  Changing the preconditioner is
    ///   allowed if the flexible option is enabled.
    void restart (const int maxIterCount);

    //! Restart with the same max iteration count as before
    void restart () { restart (maxNumIters_); }

    /// \brief "Back out" iterations following numIters 
    ///
    /// "Back out" iterations after, but not including, numIters.
    /// This is permanent, but relatively inexpensive.  The main cost
    /// is \fn$O(m^2)\fn$ floating-point operations on small dense
    /// matrices and vectors, where m = numIters().
    ///
    /// \param numIters [in] 0 <= numIters <= getNumIters().
    void backOut (const int numIters);

    /// \brief Compute and return current solution update
    ///
    /// Compute and return current solution update.  This also updates
    /// the value returned by currentNativeResidualNorm().
    Teuchos::RCP<MV> getCurrentUpdate ();

    /// \brief The current "native" residual vector
    /// 
    /// Compute the current "native" residual vector, which is formed
    /// from the basis vectors V, the upper Hessenberg matrix H, the
    /// current projected least-squares solution y, and the initial
    /// residual vector.  Computing this does not invoke the operator
    /// or preconditioner(s).
    Teuchos::RCP<MV> currentNativeResidualVector ();

    /// \brief The current "native" residual norm
    ///
    /// This is computed from the projected least-squares problem
    /// involving the upper Hessenberg matrix.  Calling this does not
    /// invoke the operator or preconditioner(s), nor does it require
    /// computations on basis vectors.
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType 
    currentNativeResidualNorm ();

    /// \brief Forward normwise error in the Arnoldi relation
    ///
    /// Compute and return the forward normwise error in the Arnoldi
    /// relation.  For a fixed preconditioner (nonflexible GMRES),
    /// we compute
    ///
    /// \fn$\| \tilde{A} V_m - V_{m+1} \underline{H}_m \|_F\fn$,
    ///
    /// where \fn$\tilde{A}\fn$ is the preconditioned matrix.  (We use
    /// the Frobenius norm to minimize dependency on the factorization
    /// codes which would be necessary for computing the 2-norm.)  For
    /// a varying right preconditioner (Flexible GMRES), we compute
    ///
    /// \fn$\| A Z_m - V_{m+1} \underline{H}_m \|_F\fn$.
    ///
    /// In these expressions, "m" is the current iteration count, as
    /// returned by getNumIters().
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType 
    arnoldiRelationError () const;

  protected:

    /// \name Subclass-specific implementation details
    ///
    /// Subclasses, by implementing all of these methods, implement a
    /// specific (F)GMRES algorithm.
    //{@

    /// \brief Whether it's legal to call extendBasis()
    ///
    /// Whether the basis (bases, if Flexible GMRES) can be extended
    /// by computing candidate basis vectors, i.e., whether it's legal
    /// to call extendBasis().  Subclasses' implementations of
    /// extendBasis() must throw GmresCantExtendBasis if no more
    /// candidate basis vectors can be computed.
    virtual bool canExtendBasis() const = 0;

    /// \brief Extend the basis by 1 (or more) vector(s)
    ///
    /// Extend the basis (bases, if Flexible GMRES) by 1 (or more)
    /// vector(s), without orthogonalizing it/them.  Standard GMRES
    /// only adds one basis vector at a time; CA-GMRES may add more
    /// than one at a time.
    ///
    /// \param V_cur [out] New basis vector(s).  The number of
    ///   column(s) gives the number of basis vector(s) added.  This
    ///   may be a view of member data, rather than freshly allocated
    ///   storage.
    ///
    /// \param Z_cur [out] If running Flexible GMRES, the
    ///   preconditioned basis vector(s); else, Teuchos::null.  The
    ///   number of column(s) gives the number of basis vector(s)
    ///   added.  This may be a view of member data, rather than
    ///   freshly allocated storage.
    ///
    /// \warning Don't call this if canExtendBasis() returns false.
    ///
    /// \note This method is non-const so that implementations may use
    ///   previously allocated space for the basis vectors.  (They are
    ///   not _required_ to use previously allocated space, so this
    ///   method may, if it wishes, allocate a new multivector to
    ///   return.)
    virtual void 
    extendBasis (Teuchos::RCP<MV>& V_cur, 
		 Teuchos::RCP<MV>& Z_cur) = 0;

    /// \brief Orthogonalize the candidate basis/es
    ///
    /// Flexible CA-GMRES must orthogonalize the Z basis along with
    /// the V basis, and it must keep both sets of orthogonalization
    /// coefficients.  Flexible standard GMRES, in contrast, only
    /// needs to orthogonalize the new V basis vector; it doesn't need
    /// to orthogonalize the Z basis vector.  This is why we let the
    /// subclass implement this functionality, rather than calling the
    /// OrthoManager's projectAndNormalize() method directly; the
    /// subclass has to decide whether it needs to orthogonalize Z_cur
    /// as well as V_cur.
    ///
    /// Subclasses may want to save the rank (if applicable) from
    /// ortho_->projectAndNormalize() somewhere.  In the case of
    /// Flexible CA-GMRES, there are _two_ ranks -- one for the V
    /// basis, and the other for the Z basis.  This is why this method
    /// doesn't return the rank, as ortho_->projectAndNormalize()
    /// does.
    ///
    /// C_V and B_V are the "C" (projection) resp. "B" (normalization)
    /// block coefficients from ortho_->projectAndNormalize() on the
    /// new V basis vectors.  If the new Z basis vectors need to be
    /// orthogonalized as well, then C_Z and B_Z are the "C" resp. "B"
    /// block coefficients from ortho_->projectAndNormalize() on the
    /// new Z basis vectors.  Subclasses' implementations of
    /// orthogonalize() are responsible for allocating space for C_V
    /// and B_V, and C_Z and B_Z if necessary (otherwise they are set
    /// to Teuchos::null, as the RCPs are passed by reference).
    ///
    /// \param V_cur [in/out] Candidate basis vector(s) for V basis
    /// \param V_prv [in] Previously orthogonalized V basis vectors
    /// \param C_V [out] Projection coefficients (V_prv^* V_cur) of
    ///   candidate V basis vector(s)
    /// \param B_V [out] Normalization coefficients (after projection)
    ///   of candidate V basis vector(s)
    ///
    /// \param Z_cur [in/out] Candidate basis vector(s) for Z basis,
    ///   or null if there is no Z basis
    /// \param Z_prv [in] Previously orthogonalized Z basis vectors,
    ///   or null if there is no Z basis
    /// \param C_Z [out] Projection coefficients (Z_prv^* Z_cur) of
    ///   candidate Z basis vector(s), or null if there is no Z basis
    /// \param B_Z [out] Normalization coefficients (after projection)
    ///   of candidate Z basis vector(s), or null if there is no Z
    ///   basis
    virtual void
    orthogonalize (const Teuchos::RCP<MV>& V_cur,
		   const Teuchos::RCP<const MV>& V_prv,
		   Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> >& C_V,
		   Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> >& B_V,
		   const Teuchos::RCP<MV>& Z_cur,
		   const Teuchos::RCP<const MV>& Z_prv,
		   Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> >& C_Z,
		   Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> >& B_Z) = 0;
    
    /// \brief Update the upper Hessenberg matrix (H_)
    ///
    /// \param C_V [in] Projection coefficients (V_prv^* V_cur) of
    ///   candidate V basis vector(s)
    /// \param B_V [in] Normalization coefficients (after projection)
    ///   of candidate V basis vector(s)
    /// \param C_Z [in] Projection coefficients (Z_prv^* Z_cur) of
    ///   candidate Z basis vector(s), or null if there is no Z basis
    /// \param B_Z [in] Normalization coefficients (after projection)
    ///   of candidate Z basis vector(s), or null if there is no Z
    ///   basis
    ///
    /// \note For Flexible GMRES, it may be desirable to compute or
    /// update a rank-revealing decomposition of the upper square
    /// submatrix of H_.  This is because FGMRES only promises
    /// convergence in exact arithmetic if this submatrix is full rank
    /// _and_ the computed residual norm is zero.  Any decomposition
    /// of H_ should _not_ be done in place; this is because \c
    /// updateProjectedLeastSquaresProblem() depends on H_ being
    /// intact.  That method does not itself modify H_, so
    /// implementations of \c updateUpperHessenbergMatrix() can also
    /// rely on H_ being intact.
    ///
    /// \note For an algorithm for updating a rank-revealing
    /// decomposition, see e.g., G. W. Stewart, "Updating a
    /// rank-revealing ULV decomposition", SIAM J. Matrix Anal. &
    /// Appl., Volume 14, Issue 2, pp. 494-499 (April 1993).
    virtual void 
    updateUpperHessenbergMatrix (const Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> >& C_V,
				 const Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> >& B_V,
				 const Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> >& C_Z,
				 const Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> >& B_Z) = 0;

    //! Whether the subclass "accepted" the candidate basis
    virtual bool acceptedCandidateBasis() const = 0;

    /// \brief Whether advance() should retry advancing the basis
    ///
    /// In advance(), if the subclass rejects the candidate basis
    /// vectors (acceptedCandidateBasis() returns false), then we have
    /// to decide whether to retry the current iteration.  Retry is
    /// done by calling advance() recursively.  Subclasses should set
    /// state before advance() calls acceptedCandidateBasis(), to
    /// ensure that this recursion terminates either with success, or
    /// by throwing GmresRejectsCandidateBasis.
    ///
    /// Default behavior is never to retry.  This always results in
    /// finite termination of advance().
    virtual bool shouldRetryAdvance() { return false; }

    //@}

    /// \name Hooks
    ///
    /// Hooks allow subclasses to preface and postface restarting,
    /// and/or the orthogonalization of new basis vector(s), with
    /// additional features.  Typical uses would be implementing
    /// Recycling GMRES (GCRO-DR), or otherwise ensuring orthogonality
    /// of the computed Krylov basis/es to another given basis.
    /// GmresBase provides trivial default implementations of all the
    /// hooks.
    ///
    /// \note In terms of aspect-oriented programming, restart() and
    ///   orthogonalize() are "join points," and the hooks are
    ///   "advice."  Unfortunately we have to implement the
    ///   "pointcuts" by hand, since otherwise we would need C++
    ///   syntax extensions.
    ///
    /// \note If you know something about Emacs Lisp, you would say
    ///   these hooks "advice" the restart() and orthogonalize()
    ///   methods.
    ///
    /// \warning Subclasses are responsible for invoking (or not
    /// invoking) their parent class' hook implementation, and for
    /// choosing when to do so.  Each child class need only invoke its
    /// immediate parent's hook; recursion will take care of the rest.
    /// However, if your class inherits from multiple classes that all
    /// inherit from a subclass of GmresBase ("diamonds" in the class
    /// hierarchy), you'll need to resolve the hooks manually.
    /// (Otherwise, recursion will invoke the same hooks multiple
    /// times.)
    //@{ 

    /// \brief Hook for subclass features before restarting.
    ///
    /// This method is called by restart(), before restart() does
    /// anything.  The default implementation does nothing.  This
    /// would be the place to add things like Krylov subspace
    /// recycling (the part where the recycling subspace is computed
    /// from the existing basis).
    ///
    /// \note Implementations of this method should not update the
    /// linear problem; this happens immediately following the
    /// invocation of preRestartHook().
    virtual void preRestartHook (const int maxIterCount) {}

    /// \brief Hook for subclass features after restarting.
    ///
    /// This method is called by restart(), after restart() finishes
    /// everything else.  The default implementation does nothing.
    /// This would be the place where subclasses may implement things
    /// like orthogonalizing the first basis vector (after restart)
    /// against the recycled subspace.
    ///
    /// \note The initial residual has already been set by restart()
    /// before this method is called.
    virtual void postRestartHook() {}

    /// \brief Hook for before orthogonalize()
    ///
    /// Hook for before calling orthogonalize() on the given
    /// newly generated basis vectors V_cur (of the V basis) and Z_cur
    /// (of the Z basis).  The default implementation does nothing.
    /// 
    /// \note Subclasses that perform subspace recycling could use
    /// implement orthogonalization against the recycled subspace
    /// either here or in \c postOrthogonalizeHook().
    virtual void 
    preOrthogonalizeHook (const Teuchos::RCP<MV>& V_cur, 
			  const Teuchos::RCP<MV>& Z_cur) {}

    /// \brief Hook for after orthogonalize()
    ///
    /// Hook for after calling orthogonalize() on the given
    /// newly generated basis vectors V_cur (of the V basis) and Z_cur
    /// (of the Z basis).  The default implementation does nothing.
    ///
    /// \note Subclasses that perform subspace recycling could use
    /// implement orthogonalization against the recycled subspace
    /// either here or in \c postOrthogonalizeHook().
    void
    postOrthogonalizeHook (const Teuchos::RCP<MV>& V_cur, 
			   const Teuchos::RCP<mat_type>& C_V,
			   const Teuchos::RCP<mat_type>& B_V,
			   const Teuchos::RCP<MV>& Z_cur,
			   const Teuchos::RCP<mat_type>& C_Z,
			   const Teuchos::RCP<mat_type>& B_Z) {}
    //@}

  private:

    //! \name Helper methods
    //@{

    /// Make a deep copy of the (left-preconditioned, if applicable)
    /// initial residual vector r0 from the linear problem instance.
    /// Compute the initial residual norm.  Allocate the V_ basis if
    /// necessary or if reallocateBasis=true.  Copy r0 into the first
    /// column of V_ and scale by the initial residual norm; that
    /// forms the first basis vector.
    ///
    /// \param maxIterCount [in] Maximum number of iterations that
    ///   GMRES may execute before restarting.  maxIterCount+1 basis
    ///   vector(s) will be allocated if necessary.
    /// \param reallocateBasis [in] If true, (re)allocate the V_
    ///   basis.  The V_ basis will be (re)allocated regardless if V_
    ///   is null or has a number of columns other than
    ///   maxIterCount+1.
    ///
    /// \warning This method should only be called by the constructor
    ///   and by restart().  It assumes that the LinearProblem
    ///   instance is set and that its setLSIndex() method has been
    ///   called with an index of exactly one right-hand side for
    ///   which to solve.
    void 
    setInitialResidual(const int maxIterCount,
		       const bool reallocateBasis = false);

    /// \brief Update the projected least-squares problem
    ///
    /// Update the QR factorization of the (CA-)GMRES upper Hessenberg
    /// matrix and its right-hand side in the projected least-squares
    /// problem, and solve the projected least-squares problem for the
    /// solution update coefficients y.
    ///
    /// \return "Native" residual norm
    ///
    /// \note Standard GMRES only updates one column at a time of the
    /// upper Hessenberg matrix; we allow updating multiple columns in
    /// order to support CA-GMRES, which computes several new columns of
    /// the upper Hessenberg matrix at a time.
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType 
    updateProjectedLeastSquaresProblem ();

    /// Views of previously orthogonalized basis vectors
    ///
    /// \param V_prv [out] View of previously orthogonalized V basis
    ///   vectors (a view into V_)
    ///
    /// \param Z_prv [out] If running Flexible GMRES, view of
    ///   previously orthogonalized Z basis vectors (a view into Z_).
    ///   Else, Teuchos::null.
    ///
    /// \warning Only call this after setInitialResidual() has been
    ///   called (that method is called in the constructor, and in
    ///   restart()).
    void 
    previousVectors (Teuchos::RCP<const MV>& V_prv,
		     Teuchos::RCP<const MV>& Z_prv) const;

    /// \brief Initial residual vector
    ///
    /// Initial residual vector (left preconditioned, if there is a
    /// left preconditioner) from the linear problem to solve.  This
    /// is copied deeply from the input linear problem.
    ///
    /// \note We keep this here so that it's easy to generate new
    ///   basis vectors, since the V_ basis comes from the same space
    ///   as the initial residual vector.
    ///
    /// \warning Please, please don't cast to nonconst MV and change
    ///   this.  You have been warned.  Make a deep copy if you want
    ///   a mutable vector.
    Teuchos::RCP<const MV> initResVec() const {
      return initResVec_;
    }

    //@}

  protected:
    //! \name Member data, inherited by subclasses
    //@{

    //! The linear problem to solve
    Teuchos::RCP<LinearProblem<Scalar, MV, OP> > lp_;

    //! Orthogonalization manager
    Teuchos::RCP<const OrthoManager<Scalar, MV> > ortho_;

    /// \brief "Native" residual vector
    ///
    /// "Native" means computed using exact-arithmetic invariants of
    /// the iterative method.  In this case, if m is the current
    /// number of iterations and \fn$r_0 = A x_0 - b\fn$,
    ///
    /// \fn$A x_k - b = A (x_0 + Z(1:m) y(1:m)) - b 
    ///               = r_0 + A Z(1:m) y(1:m) 
    ///               = r_0 + V(1:m+1) H(1:m+1, 1:m) y(1:m)\fn$.
    ///
    /// Accuracy of the above formula depends only on the Arnoldi
    /// relation, which in turn depends mainly on the residual error
    /// of the orthogonalization being small.  This is true for just
    /// about any orthogonalization method, so the above computation
    /// should hold just about always.
    ///
    /// Storage for the native residual vector is allocated only on
    /// demand.  The storage is cached and reallocated only when the
    /// LinearProblem changes.  (This ensures that the dimensions and
    /// data distribution are correct, i.e., are the same as the
    /// initial residual vector in the linear problem).
    Teuchos::RCP<MV> nativeResVec_;

    /// \brief Initial residual vector.
    ///
    /// This is the left-preconditioned residual vector if a left
    /// preconditioner has been provided, otherwise it is the
    /// unpreconditioned residual vector.
    Teuchos::RCP<MV> initResVec_;
    
    /// \brief Current solution update
    ///
    /// The current approximate solution for Flexible GMRES is 
    /// x0 + Z_[0:numIters-1]y, and for (nonflexible) GMRES is
    /// x0 + V_[0:numIters-1]y.
    ///
    /// \note Invalidated when linear problem changes, etc.
    Teuchos::RCP<MV> xUpdate_;

    /// \brief First set of (orthogonalized) Krylov basis vectors
    ///
    /// \note These vectors are always used, whether we are performing
    /// Flexible GMRES (left preconditing with a possibly changing
    /// preconditioner, so we keep a second set of basis vectors Z_)
    /// or standard GMRES.
    Teuchos::RCP<MV> V_;

    /// \brief Second set of (orthogonalized) Krylov basis vectors
    ///
    /// \note These vectors are only used if we are performing
    /// Flexible GMRES.  In that case, the solution update
    /// coefficients are computed for this basis.
    Teuchos::RCP<MV> Z_;

    /// \brief The Arnoldi/GMRES upper Hessenberg matrix
    ///
    /// H_[0:curNumIters_, 0:curNumIters_-1] (inclusive index ranges,
    /// not exclusive like SciPy's) is always valid and is the upper
    /// Hessenberg matrix for the first curNumIters_ iterations of
    /// Arnoldi/GMRES.  
    ///
    /// \note We do not overwrite H_ when computing its QR
    ///   factorization; new H_ data is copied into R_ and updated in
    ///   place there.
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int, Scalar> > H_;

    //! The R factor in the QR factorization of H_
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int, Scalar> > R_;

    //! Current solution of the projected least-squares problem
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int, Scalar> > y_;

    /// \brief Current RHS of the projected least-squares problem
    ///
    /// The current right-hand side of the projected least-squares
    /// problem \fn$\min_y \|\underline{H} y - \beta e_1\|_2\fn$.  z_
    /// starts out as \fn$\beta e_1\fn$ (where \fn$\beta\fn$ is the
    /// initial residual norm).  It is updated progressively along
    /// with the QR factorization of H_.
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int, Scalar> > z_;

    Teuchos::Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> theCosines_;
    Teuchos::Array<Scalar> theSines_;

    /// \brief The initial residual norm
    ///
    /// GMRES makes use of the initial residual norm for solving the
    /// projected least-squares problem for the solution update
    /// coefficients.  In that case, it is usually called
    /// \fn$\beta\fn$.  For left-preconditioned GMRES, this is the
    /// preconditioned initial residual norm, else it's the
    /// unpreconditioned version.
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType initialResidualNorm_;

    /// Last column of H_ for which the QR factorization (implicitly
    /// stored in theCosines_, theSines_, and R_) has been computed.
    /// Should be initialized to -1 (which means the 0th column of H_
    /// will be the first to be updated).
    int lastUpdatedCol_;

    /// \brief Current number of completed iterations.
    ///
    /// This also counts the number of orthogonalized basis vectors in
    /// V_.  Columns [0, curNumIters_] of V_ are always orthogonalized
    /// basis vectors.
    int curNumIters_;

    /// \brief Maximum number of iterations allowed
    ///
    /// This affects the amount of storage consumed by V_ and Z_.  Our
    /// Arnoldi/GMRES implementation requires that all the basis
    /// vectors in a set occupy a contiguous block of storage.  (This
    /// differs from some other Belos implementations of iterative
    /// methods, where each basis vector (or block of vectors) is
    /// allocated separately on demand.)  
    int maxNumIters_;

    /// \brief Whether we are running Flexible GMRES
    ///
    /// Flexible GMRES (FGMRES) is a variant of standard GMRES that
    /// allows the preconditioner to change in every iteration.  It
    /// only works for right preconditioning.  FGMRES requires keeping
    /// two sets of Krylov basis vectors: one for the Krylov subspace
    ///
    /// \fn$\text{span}\{r_0, A M^{-1} r_0, \dots, (A M^{-1})^k r_0\}\fn$ 
    ///
    /// where \fn$r_0\fn$ is the unpreconditioned residual, and one
    /// for the Krylov subspace
    ///
    /// \fn$\text{span}\{M^{-1} r_0, M^{-1} A M^{-1} r_0, \dots, M^{-1} (A M^{-1})^{k-1} r_0\}\fn$.
    ///
    /// We store the basis for the first subspace in V_, and the basis
    /// for the second subspace in Z_.  If we are not running FGMRES,
    /// we let Z_ be Teuchos::null.  FGMRES reduces to standard GMRES
    /// in the case of nopreconditioning, and then we also let Z_ be
    /// null.
    ///
    /// \note In the original Flexible GMRES paper, Saad suggested
    /// implementing FGMRES and GMRES in a single routine, with a
    /// runtime switch to decide whether to keep the Z basis and
    /// whether update the solution with the Q or Z basis.  This is in
    /// fact what we do.  The differences between GMRES and FGMRES are
    /// small enough that we can benefit from the reduced code
    /// duplication.  We have smart algorithms people, but not many of
    /// them (and they are very busy), so it's a more efficient use of
    /// development effort to have less code that is slightly more
    /// complicated.
    bool flexible_;

    //@}
  };

  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  template<class Scalar, class MV, class OP>
  void
  GmresBase<Scalar, MV, OP>::
  previousVectors (Teuchos::RCP<const MV>& V_prv,
		   Teuchos::RCP<const MV>& Z_prv) const
  {
    using Teuchos::Range1D;
    
    // The number of iterations (returned by getNumIters()) does not
    // include the initial basis vector, which is the first column of
    // V_.
    const int m = getNumIters(); 
    V_prv = MVT::CloneView (*V_, Range1D(0,m));
    if (flexible_)
      Z_prv = MVT::CloneView (*Z_, Range1D(0,m));
    else
      Z_prv = Teuchos::null;
  }

  template<class Scalar, class MV, class OP>
  void
  GmresBase<Scalar, MV, OP>::advance () 
  {
    using Teuchos::is_null;
    using Teuchos::null;
    using Teuchos::Range1D;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::tuple;
    typedef Teuchos::SerialDenseMatrix<int,Scalar> mat_type;

    TEST_FOR_EXCEPTION( !canExtendBasis(), GmresCantExtendBasis,
			"GMRES (iteration " << getNumIters() << ") cannot "
			"add any more basis vectors." );

    RCP<const MV> V_prv, Z_prv;
    RCP<MV> V_cur, Z_cur;
    previousVectors (V_prv, Z_prv);
    extendBasis (V_cur, Z_cur);
    preOrthogonalizeHook (V_cur, Z_cur);

    // Flexible CA-GMRES must orthogonalize the Z basis along with the
    // V basis, and it must keep both sets of orthogonalization
    // coefficients.  Flexible standard GMRES only needs to
    // orthogonalize the new V basis vector; it doesn't need to
    // orthogonalize the Z basis vector.  This is why we let the
    // subclass implement this functionality, rather than calling the
    // OrthoManager's projectAndNormalize() method directly; the
    // subclass has to decide whether it needs to orthogonalize Z_cur
    // as well as V_cur.
    //
    // Subclasses may want to save the rank (if applicable) from
    // ortho_->projectAndNormalize() somewhere.  In the case of
    // Flexible CA-GMRES, there are _two_ ranks -- one for the V
    // basis, and the other for the Z basis.
    //
    // C_V and B_V are the "C" (projection) resp. "B" (normalization)
    // block coefficients from ortho_->projectAndNormalize() on the
    // new V basis vectors.  If the new Z basis vectors need to be
    // orthogonalized as well, then C_Z and B_Z are the "C" resp. "B"
    // block coefficients from ortho_->projectAndNormalize() on the
    // new Z basis vectors.  Subclasses' implementations of
    // orthogonalize() are responsible for allocating space for
    // C_V and B_V, and C_Z and B_Z if necessary (otherwise they are
    // set to Teuchos::null, as the RCPs are passed by reference).
    RCP<mat_type> C_V, B_V, C_Z, B_Z;
    orthogonalize (V_cur, V_prv, C_V, B_V, Z_cur, Z_prv, C_Z, B_Z);
    postOrthogonalizeHook (V_cur, C_V, B_V, Z_cur, C_Z, B_Z);
    // Implemented by subclasses.  C_Z and B_Z may not be used (in
    // which case orthogonalize() may set them to Teuchos::null).
    updateUpperHessenbergMatrix (C_V, B_V, C_Z, B_Z);

    // The subclass decides whether or not to accept the candidate
    // basis.
    //
    // We update the upper Hessenberg matrix before deciding
    // acceptance, because the acceptance critera might depend on the
    // entire upper Hessenberg matrix (e.g., (an estimate of) its
    // condition number), not just the condition number or rank of the
    // candidate basis vector(s).
    //
    // Subclasses' implementations of orthogonalize() and/or
    // updateUpperHessenbergMatrix() might wish to set state to
    // indicate acceptance or rejection of the candidate basis
    // vector(s).  They may also call advance() recursively (e.g.,
    // setting the "last rejected rank" so that the next call of
    // advance() computes fewer candidate basis vectors).
    if (! acceptedCandidateBasis())
      rejectCandidateBasis ();
    else
      acceptCandidateBasis ();

    // We don't update the projected least-squares problem here.
    // This is done lazily, whenever the "native" residual norm
    // or the current solution update are requested.
  }


  template<class Scalar, class MV, class OP>
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType
  GmresBase<Scalar, MV, OP>::arnoldiRelationError () const
  {
    using Teuchos::Range1D;
    using Teuchos::RCP;
    const Scalar one = STS::one();
    const int m = getNumIters();

    if (m == 0) // No iterations completed means no error (yet)
      return magnitude_type(0);
    RCP<const MV> V_mp1 = MVT::CloneView (*V_, Range1D(0,m));
    const mat_type H_m (Teuchos::View, *H_, m+1, m);
    std::vector<magnitude_type> norms (m);
    if (flexible_)
      {
	RCP<const MV> Z_view = MVT::CloneView (*Z_, Range1D(0,m-1));
	// A*Z is in the same space as V, so we clone from V
	RCP<MV> AZ = MVT::Clone (*V_, m);
	// Compute AZ := A*Z_m
	lp_->applyOp (Z_view, AZ);
	// Compute AZ := AZ - V_{m+1} \underline{H}_m
	MVT::MvTimesMatAddMv (-one, V_mp1, H_m, one, AZ);
	MVT::MvNorm (AZ, norms);
      }
    else
      {
	RCP<const MV> V_m = MVT::CloneView (*V_, Range1D(0,m-1));
	RCP<MV> AV = MVT::Clone (*V_, m);
	// Compute AV := A*V_m.  Here, "A" means the preconditioned
	// operator.
	lp_->apply (V_m, AV);
	// Compute AV := A V_m - V_{m+1} \underline{H}_m
	MVT::MvTimesMatAddMv (-one, V_mp1, H_m, one, AV);
	MVT::MvNorm (AV, norms);
      }
    // Compute and return the Frobenius norm of the above (either AZ
    // or AV, depending on whether we are performing Flexible GMRES).
    //
    // FIXME (mfh 16 Feb 2011) This computation might overflow for
    // unjustifiable reasons.  I should really include intermediate
    // scaling.
    magnitude_type result (0);
    for (int k = 0; k < m; ++k)
      result += norms[k];
    return STM::squareroot(result);
  }


  template<class Scalar, class MV, class OP>
  GmresBase<Scalar, MV, OP>::
  GmresBase (const Teuchos::RCP<LinearProblem<Scalar, MV, OP> >& problem,
	     const Teuchos::RCP<const OrthoManager<Scalar, MV> >& ortho,
	     const int maxIterCount,
	     const bool flexible) :
    lp_ (problem), 
    ortho_ (ortho),
    lastUpdatedCol_ (-1), // column updates start with zero
    curNumIters_ (0),
    maxNumIters_ (maxIterCount),
    flexible_ (flexible && ! lp_.is_null() && ! lp_->getRightPrec().is_null())
  {
    using Teuchos::is_null;
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::rcp;

    // Fragments of error messages for use in sanity checks.
    const char prefix[] = "Belos::GmresBase constructor: ";
    const char postfix[] = "  Please call setProblem() with non-null data "
      "before passing the LinearProblem to the GmresBase constructor.";
    const char postfix2[] = "  That may mean that the LinearProblem's "
      "setLSIndex() method has not yet been called, even though setProblem() "
      "has been called.  After calling setProblem(), you should call "
      "setLSIndex() with the ind{ex,ices} of the right-hand side(s) to solve.";

    // Sanity checks on the linear problem.
    //
    // First, make sure A, X, and B are all non-null.
    TEST_FOR_EXCEPTION(lp_.is_null(), std::invalid_argument,
		       prefix << "The given LinearProblem is null.");
    TEST_FOR_EXCEPTION(! lp_->isProblemSet(), std::invalid_argument,
    		       prefix << "The given LinearProblem has not yet been set."
		       << postfix);
    TEST_FOR_EXCEPTION(lp_->getLHS().is_null(), std::invalid_argument,
    		       prefix << "The given LinearProblem has null initial guess"
		       " (getLHS()) for all right-hand sides." << postfix);
    TEST_FOR_EXCEPTION(lp_->getRHS().is_null(), std::invalid_argument,
    		       prefix << "The given LinearProblem's right-hand side(s) "
		       "are all null (getRHS().is_null() == true)." << postfix);
    TEST_FOR_EXCEPTION(lp_->getOperator().is_null(), std::invalid_argument,
    		       prefix << "The given LinearProblem's operator (the "
		       "matrix A) is null." << postfix);
    // Next, make sure that setLSIndex() has been called on the linear
    // problem instance, so that the "current" right-hand side and
    // "current" initial guess have been set.
    TEST_FOR_EXCEPTION(lp_->getCurrLHSVec().is_null(), 
		       std::invalid_argument,
    		       prefix << "Although the given LinearProblem has non-null "
		       "initial guess (getLHS()) for all right-hand sides, the "
		       "current initial guess (getCurrLHSVec()) is null."
		       << postfix2);
    TEST_FOR_EXCEPTION(lp_->getCurrRHSVec().is_null(), 
		       std::invalid_argument,
    		       prefix << "Although the given LinearProblem has non-null "
		       "initial guess (getLHS()) for all right-hand sides, the "
		       "current right-hand side to solve (getCurrRHSVec()) is "
		       "null." << postfix2);
    // Get the initial guess x0, and allocate the solution update
    // (xUpdate_) from the same vector space as x0.
    RCP<const MV> x0 = lp_->getCurrLHSVec();
    {
      const int numLHS = MVT::GetNumberVecs (*x0);
      TEST_FOR_EXCEPTION(numLHS != 1, std::invalid_argument,
			 "Our GMRES implementation only works for single-"
			 "vector problems, but the supplied initial guess has "
			 << numLHS << " columns.");
      const int numRHS = MVT::GetNumberVecs (*(lp_->getCurrRHSVec()));
      TEST_FOR_EXCEPTION(numRHS != 1, std::invalid_argument,
			 "Our GMRES implementation only works for single-"
			 "vector problems, but the current right-hand side has "
			 << numRHS << " columns.");
    }
    xUpdate_ = MVT::Clone (*x0, 1);
    MVT::MvInit (*xUpdate_, STS::zero());

    // Get the (left preconditioned, if applicable) residual vector,
    // and make a deep copy of it in initResVec_.  Allocate the V_
    // basis vectors, compute the initial residual norm, and compute
    // the first basis vector (first column of V_).
    setInitialResidual (maxIterCount);

    // The Z_ vectors, if we need them, are always in the same vector
    // space as the initial solution guess (x0).  Even if the caller
    // asks for the flexible option, we don't allocate space for Z_
    // unless the caller specifies a right preconditioner.
    Z_ = (flexible && ! lp_->getRightPrec().is_null()) ? 
      MVT::Clone(*x0, maxIterCount+1) : null;

    // These (small dense) matrices and vectors encode the projected
    // least-squares problem.  z_ was already allocated and
    // initialized in setInitialResidual() above.
    H_ = rcp (new mat_type (maxIterCount+1, maxIterCount));
    R_ = rcp (new mat_type (maxIterCount+1, maxIterCount));
    y_ = rcp (new mat_type (maxIterCount, 1));

    // These cosines and sines encode the Q factor in the QR
    // factorization of the upper Hessenberg matrix H_.  We compute
    // this by computing H_ into R_ and operating on R_ in place; H_
    // itself is left alone.
    theCosines_.resize (maxIterCount);
    theSines_.resize (maxIterCount);
  }

  template<class Scalar, class MV, class OP>
  Teuchos::RCP<MV> 
  GmresBase<Scalar, MV, OP>::currentNativeResidualVector ()
  {
    const Scalar zero = STS::zero();
    const Scalar one = STS::one();
    using Teuchos::is_null;
    using Teuchos::Range1D;
    using Teuchos::RCP;

    const char prefix[] = "Belos::GmresBase::currentNativeResidualVector: ";
    TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*initResVec_) != 1,
		       std::logic_error,
		       prefix << "Initial residual vector (initResVec_) has " 
		       << MVT::GetNumberVecs(*initResVec_)
		       << " columns, but should only have 1 column.  "
		       "This should never happen.");
    // If left preconditioning is used, then the "native" residual
    // vector is in the same vector space as the left-preconditioned
    // initial residual vector.  Otherwise, it is in the same space as
    // the right-hand side b and the unpreconditioned initial residual
    // vector.
    //
    // Assign initial residual vector to nativeResVec; we will update
    // it below if the number of iterations is > 0.  ("Native"
    // residual vector is the same as the initial residual vector when
    // the number of iterations is 0.)
    if (nativeResVec_.is_null() || MVT::GetNumberVecs (*nativeResVec_) != 1)
      // Allocate nativeResVec if necessary.
      nativeResVec_ = MVT::CloneCopy (*initResVec_, Range1D(0,0));
    else
      MVT::SetBlock (*initResVec_, Range1D(0,0), *nativeResVec_);

    if (curNumIters_ > 0)
      {
	TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*V_) < curNumIters_+1, 
			   std::logic_error,
			   "Only " << MVT::GetNumberVecs(*V_) << " basis vectors "
			   "were given, but curNumIters+1=" << (curNumIters_+1) 
			   << " of them are required.  "
			   "This likely indicates a bug in Belos.");
	TEST_FOR_EXCEPTION(H_->numRows() < curNumIters_+1, 
			   std::logic_error,
			   "H only has " << H_->numRows() << " rows, but " 
			   << (curNumIters_+1) << " rows are required.  "
			   "This likely indicates a bug in Belos.");
	TEST_FOR_EXCEPTION(H_->numCols() < curNumIters_, 
			   std::logic_error,
			   "H only has " << H_->numCols() << " columns, but "
			   << (curNumIters_+1) << " columns are required.  "
			   "This likely indicates a bug in Belos.");
	// TEST_FOR_EXCEPTION(y_->length() < curNumIters_, std::logic_error,
	// 		   "y only has " << y_->length() << " entries, but "
	// 		   curNumIters_ << " entries are required.  "
	// 		   "This likely indicates a bug in Belos.");
	RCP<const MV> V_view = MVT::CloneView (*V_, Range1D(0, curNumIters_));
	const mat_type H_view (Teuchos::View, *H_, curNumIters_+1, curNumIters_);
	const mat_type y_view (Teuchos::View, *y_, curNumIters_, 1);
	mat_type H_times_y (curNumIters_+1, 1);
	{
	  const int err = 
	    H_times_y.multiply (Teuchos::NO_TRANS, Teuchos::NO_TRANS, 
				one, H_view, y_view, zero);
	  TEST_FOR_EXCEPTION(err != 0, std::logic_error,
			     "In GMRES, when computing the current native resi"
			     "dual vector via the Arnoldi relation, H*y failed"
			     " due to incompatible dimensions.  This is likely"
			     " a Belos bug.");
	}
	// nativeResVec_ := V * (H * y) - r0 (where we stored r0 in nativeResVec_)
	MVT::MvTimesMatAddMv (one, *V_view, H_times_y, -one, *nativeResVec_);
      }
    return nativeResVec_;
  }

  template<class Scalar, class MV, class OP>
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType
  GmresBase<Scalar, MV, OP>::currentNativeResidualNorm ()
  {
    // Update the least-squares problem if necessary.  This computes
    // the current native residual norm for "free."
    return updateProjectedLeastSquaresProblem();
  }

  template<class Scalar, class MV, class OP>
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType
  GmresBase<Scalar, MV, OP>::updateProjectedLeastSquaresProblem()
  {
    const scalar_type zero = STS::zero();
    const scalar_type one = STS::one();
    const int startCol = lastUpdatedCol_ + 1;
    const int endCol = curNumIters_ - 1;
    blas_type blas;

    TEST_FOR_EXCEPTION(startCol > endCol+1, std::logic_error,
		       "Somehow, GmresBase updated the QR factorization of the "
		       "upper Hessenberg matrix past the last updated column.  "
		       "The last updated column was " << lastUpdatedCol_ << ", "
		       "but the rightmost column to update is " << endCol << "."
		       "  This is likely a bug in GmresBase.");
    if (startCol == endCol+1) // Nothing to do
      return STS::magnitude ((*z_)(endCol+1,0));

    const int numCols = endCol - startCol + 1;
    {
      // Copy columns [startCol, endCol] (inclusive) of the current
      // upper Hessenberg matrix H, into the upper triangular matrix
      // R.  (We avoid overwriting H because having the original upper
      // Hessenberg matrix around can save us some work when computing
      // the "native" residual.)  Columns [0, startCol-1] of R are
      // already upper triangular; we don't have to do any more work
      // there.
      const int numRows = endCol + 2; // include last subdiagonal element
      const mat_type H_view (Teuchos::View, *H_, numRows, numCols, 0, startCol);
      mat_type R_view (Teuchos::View, *R_, numRows, numCols, 0, startCol);
      R_view.assign (H_view);
    }
    // Define some references to ease notation.
    mat_type& R = *R_;
    mat_type& z = *z_;
    
    // Update columns [startCol, endCol] of the QR factorization of
    // the upper Hessenberg matrix.  Use Givens rotations to maintain
    // the factorization.  We use a left-looking update procedure,
    // since Arnoldi is always adding new column(s) on the right of
    // the upper Hessenberg matrix.
    const int LDR = R.stride();
    for (int j = 0; j < startCol; ++j)
      {
	// Apply all previous Givens rotations to the new columns of
	// the upper Hessenberg matrix (as stored in R).
	//
	// FIXME (mfh 25 Dec 2010) Teuchos::BLAS wants nonconst
	// pointers for the sine and cosine arguments to ROT().  They
	// should be const, since ROT doesn't change them (???).  Fix
	// this in Teuchos.
	magnitude_type theCosine = theCosines_[j];
	scalar_type theSine = theSines_[j];
	blas.ROT (numCols, &R(j,j), LDR, &R(j+1,j), LDR, 
		  &theCosine, &theSine);
      }
    for (int j = startCol; j != endCol; ++j)
      {
	// Calculate new Givens rotation for [R(j,j); R(j+1,j)]
	magnitude_type theCosine;
	scalar_type theSine;
	blas.ROTG (&R(j,j), &R(j+1,j), &theCosine, &theSine);
	theCosines_[j] = theCosine;
	theSines_[j] = theSine;
	// Clear the subdiagonal element R(j+1,j)
	R(j+1,j) = zero;
	// Update the "trailing matrix."  The "if" check is not
	// strictly necessary, but ensures that the references to
	// column j+1 of R always makes sense (even though ROT
	// shouldn't access column j+1 if endCol-1-j==0).
	if (j < endCol)
	  blas.ROT (endCol-1-j, &R(j,j+1), LDR, &R(j+1,j+1), LDR, 
		    &theCosine, &theSine);
	// Update the right-hand side of the least-squares problem
	// with the new Givens rotation.
	blas.ROT (1, &z(j,0), 1, &z(j+1,0), 1, &theCosine, &theSine);
      }
    // Remember the rightmost updated column.
    lastUpdatedCol_ = endCol;

    // The absolute value of the last element of z gives the current
    // "native" residual norm.
    const magnitude_type nativeResNorm = STS::magnitude ((*z_)(endCol+1,0));

    // Now that we have the updated R factor of H, and the updated
    // right-hand side z, solve the least-squares problem by solving
    // the linear system Ry=z.
    // 
    // TRSM() overwrites the right-hand side with the solution, so
    // copy z into y.
    {
      mat_type y_view (Teuchos::View, *y_, endCol+1, 1);
      mat_type z_view (Teuchos::View, *z_, endCol+1, 1);
      y_view.assign (z_view);
    }
    // Solve Ry = z for y.
    blas.TRSM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
	      Teuchos::NON_UNIT_DIAG, curNumIters_, 1, 
	      one, R_->values(), LDR, y_->values(), y_->stride());
    return nativeResNorm;
  }


  template<class Scalar, class MV, class OP>
#if 0
  std::pair< Teuchos::RCP<MV>, typename Teuchos::ScalarTraits<Scalar>::magnitudeType >
#else // not 0
  Teuchos::RCP<MV>
#endif // 0
  GmresBase<Scalar, MV, OP>::getCurrentUpdate ()
  {
    using Teuchos::is_null;
    using Teuchos::Range1D;
    using Teuchos::RCP;
    const Scalar one = STS::one();
    const Scalar zero = STS::zero();

    if (is_null (xUpdate_))
      {
	// xUpdate_ comes from the Z_ space and not the V_ space in
	// the flexible variant of GMRES.  This makes a difference
	// only if there are multiple vector spaces involved, that
	// is if the preconditioner maps V-space vectors to Z-space
	// vectors, and the matrix maps Z-space vectors to V-space
	// vectors (where V and Z are different spaces).  This is
	// legitimate as long as the initial guess for the solution
	// (x0) comes from the Z space.
	if (flexible_)
	  xUpdate_ = MVT::Clone (*Z_, 1);
	else
	  xUpdate_ = MVT::Clone (*V_, 1);
      }
    if (curNumIters_ == 0)
      {
	MVT::MvInit (*xUpdate_, zero);
	return xUpdate_;
      }

    // When one iteration of GMRES is completed, the upper
    // Hessenberg matrix H is 2 by 1.  Thus, we want to update
    // columns up to and including curNumIters_ - 1: actually,
    // columns [lastUpdatedCol_+1, curNumIters_-1] (inclusive).  The
    // solution update gives you the current "native" residual norm
    // for free.  We could save and return that as well...
    //const magnitude_type nativeResNorm = updateProjectedLeastSquaresProblem();
    (void) updateProjectedLeastSquaresProblem();

    // Current iteration count does not include the initial basis vector.
    const int m = getNumIters();
    // Basis for solution update coefficients; Flexible GMRES uses the
    // preconditioned basis here.
    RCP<const MV> Z_view = flexible_ ? 
      MVT::CloneView (*Z_, Range1D(0, m-1)) : 
      MVT::CloneView (*V_, Range1D(0, m-1));
    const mat_type y_view (Teuchos::View, *y_, m, 1);
    // xUpdate_ := Z_view * y_view
    MVT::MvTimesMatAddMv (one, *Z_view, y_view, zero, *xUpdate_);
    return xUpdate_;
  }

  template<class Scalar, class MV, class OP>
  void
  GmresBase<Scalar, MV, OP>::
  setInitialResidual(const int maxIterCount,
		     const bool reallocateBasis)
  {
    using Teuchos::Range1D;
    using Teuchos::RCP;

    // Message fragment for error messages.
    const char prefix[] = "Belos::GmresBase::setInitialResidual: ";

    TEST_FOR_EXCEPTION(maxIterCount < 0, std::invalid_argument,
		       prefix << "maxIterCount = " << maxIterCount 
		       << " < 0.");
    // Do we have a left preconditioner?
    const bool haveLeftPrec = ! lp_->getLeftPrec().is_null();

    // mfh 22 Feb 2011: LinearProblem has the annoying quirk that
    // getInit(Prec)ResVec() return the residual vector for the whole
    // linear problem (X_ and B_), not for the "current" linear
    // problem (curX_ and curB_, returned by getCurrLHSVec()
    // resp. getCurrRHSVec()).  In order to get the residual vector(s)
    // for the _current_ linear problem, we have to use getLSIndex()
    // to extract the column ind(ex/ices) and take the corresponding
    // column of getInit(Prec)ResVec().  It should be exactly one
    // column, since this GMRES solver only knows how to solve for one
    // right-hand side at a time.
    RCP<const MV> R_full = 
      haveLeftPrec ? lp_->getInitPrecResVec() : lp_->getInitResVec();
    std::vector<int> inputIndices = lp_->getLSIndex();
    TEST_FOR_EXCEPTION(inputIndices.size() == 0, 
		       std::invalid_argument,
		       "The LinearProblem claims that there are zero linear "
		       "systems to solve: getLSIndex() returns an index "
		       "vector of length zero.");

    // Some of the indices returned by getLSIndex() might be -1.
    // These are for column(s) of getCurr{L,R}HSVec() that don't
    // correspond to any actual linear system(s) that we want to
    // solve.  They may be filled in by block iterative methods in
    // order to make the blocks full rank or otherwise improve
    // convergence.  We shouldn't have any of those columns here, but
    // we check just to make sure.
    std::vector<int> outputIndices;
    std::remove_copy_if (inputIndices.begin(), inputIndices.end(), 
			 std::back_inserter (outputIndices), 
			 std::bind2nd (std::equal_to<int>(), -1));

    // outputIndices should have exactly one entry, corresponding to
    // the one linear system to solve.  Our GMRES implementation can
    // only solve for one right-hand side at a time.
    if (outputIndices.size() != 1 || outputIndices[0] == -1)
      {
	std::string inputIndicesAsString;
	{
	  std::ostringstream os;
	  os << "[";
	  std::copy (inputIndices.begin(), inputIndices.end(),
		     std::ostream_iterator<int>(os, " "));
	  os << "]";
	  inputIndicesAsString = os.str();
	}
	std::string outputIndicesAsString;
	{
	  std::ostringstream os;
	  os << "[";
	  std::copy (outputIndices.begin(), outputIndices.end(),
		     std::ostream_iterator<int>(os, " "));
	  os << "]";
	  outputIndicesAsString = os.str();
	}
	std::ostringstream os;
	os << prefix << "The LinearProblem instance's getLSIndex() method "
	  "returns indices " << inputIndicesAsString << ", of which the "
	  "following are not -1: " << outputIndicesAsString << ".";
	TEST_FOR_EXCEPTION(outputIndices.size() != 1,
			   std::invalid_argument, 
			   os.str() << "  The latter list should have length "
			   "exactly 1, since our GMRES implementation only "
			   "knows how to solve for one right-hand side at a "
			   "time.");
	TEST_FOR_EXCEPTION(outputIndices[0] == -1,
			   std::invalid_argument,
			   os.str() << "  the latter list contains no "
			   "nonnegative entries, meaning that there is no "
			   "valid linear system to solve.");
      }
    // View of the "current" linear system's residual vector.
    RCP<const MV> r0 = MVT::CloneView (*R_full, outputIndices);
    // Sanity check that the residual vector has exactly 1 column.
    // We've already checked this a different way above.
    const int numResidVecs = MVT::GetNumberVecs (*r0);
    TEST_FOR_EXCEPTION(numResidVecs != 1, std::logic_error,
		       prefix << "Residual vector has " << numResidVecs 
		       << " columns, but should only have one.  "
		       "Should never get here.");
    // Save a deep copy of the initial residual vector.
    initResVec_ = MVT::CloneCopy (*r0);

    // Allocate space for the V_ basis vectors if necessary.
    //
    // If left preconditioning is used, then the V_ vectors are in the
    // same vector space as the (left-)preconditioned initial residual
    // vector.  Otherwise they are in the same space as the
    // unpreconditioned initial residual vector.  That's r0 in either
    // case.
    //
    // FIXME (mfh 22 Feb 2011) We should really check to make sure
    // that the columns of V_ (if previously allocated) are in the
    // same vector space as r0.  Belos::MultiVecTraits doesn't give us
    // a way to do that.  Just checking if the number of rows is the
    // same is not enough.  The third Boolean argument
    // (reallocateBasis) lets us add in the check externally if
    // MultiVecTraits gains that ability sometime.
    if (V_.is_null() || 
	MVT::GetNumberVecs(*V_) != maxIterCount+1 || 
	reallocateBasis)
      V_ = MVT::Clone (*initResVec_, maxIterCount+1);

    // Initial residual norm is computed with respect to the inner
    // product defined by the OrthoManager.  Since the basis
    // vectors are orthonormal (resp. unitary) with respect to that
    // inner product, this ensures that the least-squares minimization
    // is also computed with respect to that inner product.
    std::vector<magnitude_type> result (1);
    ortho_->norm (*r0, result);
    initialResidualNorm_ = result[0];

    // z_ is the right-hand side of the projected least-squares
    // problem.  If we reallocated z_, we have to reset it by setting
    // its first entry to the initial residual norm.
    if (z_.is_null() || z_->numRows() < maxIterCount+1)
      z_ = rcp (new mat_type (maxIterCount+1, 1));
    else
      (void) z_->putScalar (STS::zero());
    (*z_)(0,0) = Scalar (initialResidualNorm_);

    // Copy the initial residual vector into the first column of V_,
    // and scale the vector by its norm.
    MVT::SetBlock (*initResVec_, Range1D(0,0), *V_);
    RCP<MV> v1 = MVT::CloneViewNonConst (*V_, Range1D(0,0));
    MVT::MvScale (*v1, Scalar(1)/initialResidualNorm_);
  }

  template<class Scalar, class MV, class OP>
  void 
  GmresBase<Scalar, MV, OP>::backOut (const int numIters)
  {
    const Scalar zero = STS::zero();

    TEST_FOR_EXCEPTION(numIters < 0, std::invalid_argument,
		       "The GMRES iteration count cannot be less than "
		       "zero, but you specified numIters = " << numIters);
    TEST_FOR_EXCEPTION(numIters > getNumIters(), std::invalid_argument,
		       "The GMRES iteration count cannot be reset to more than "
		       "its current iteration count " << getNumIters()
		       << ", but you specified numIters = " << numIters << ".");
    // Reset the R factor in the QR factorization of H_, and the
    // right-hand side z_ of the projected least-squares problem.  We
    // will "replay" the first numIters steps of the QR factorization
    // below.
    //
    // The integer flag return value of SerialDenseMatrix's
    // putScalar() is not informative; it's always zero.
    (void) R_->putScalar(zero); // zero out below-diagonal entries
    *R_ = *H_; // deep copy
    (void) y_->putScalar(zero);
    (void) z_->putScalar(zero);
    (*z_)(0,0) = Scalar (initialResidualNorm_);

    // Replay the first numIters-1 Givens rotations.
    lastUpdatedCol_ = -1;
    curNumIters_ = numIters;
    (void) updateProjectedLeastSquaresProblem ();
  }

  template<class Scalar, class MV, class OP>
  void 
  GmresBase<Scalar, MV, OP>::
  restart (const int maxIterCount)
  {
    using Teuchos::RCP;
    const Scalar zero = STS::zero();

    // This would be the place where subclasses may implement things
    // like recycling basis vectors for the next restart cycle.  Note
    // that the linear problem hasn't been updated yet; this happens
    // below.
    preRestartHook (maxIterCount);

    // Make sure that the LinearProblem object gets the current
    // approximate solution (which becomes the initial guess for the
    // upcoming new restart cycle), so that the "exact" residuals are
    // correct.
    (void) getCurrentUpdate (); // results stored in xUpdate_
    (void) lp_->updateSolution (xUpdate_, true, STS::one());

    // Check whether there is a right preconditioner, since (at least
    // in theory) the caller could have added a right preconditioner
    // after a restart cycle.  (Is that legitimate?  Adding or
    // changing a right preconditioner between restarts wouldn't
    // change the initial residual vector, and therefore wouldn't
    // change the V_ basis.  The V_ basis is all that matters for the
    // restart.)
    const bool needToAllocateZ = flexible_ && ! lp_->getRightPrec().is_null();

    // This (re)allocates the V_ basis if necessary.  It also
    // (re)allocates the projected least-squares problem's right-hand
    // side z_ if necessary, and initializes it regardless.
    setInitialResidual (maxIterCount);

    // If the new max iteration count has changed, reallocate space
    // for the Z_ basis vectors and the projected least-squares
    // problem.  setInitialResidual() already does that for the V_
    // basis.
    if (maxNumIters_ != maxIterCount)
      {
	maxNumIters_ = maxIterCount;
	// Initial guess.
	RCP<MV> x0 = lp_->getCurrLHSVec();
	Z_ = needToAllocateZ ? MVT::Clone (*x0, maxIterCount+1) : null;
	H_ = rcp (new mat_type (maxIterCount+1, maxIterCount));
	R_ = rcp (new mat_type (maxIterCount+1, maxIterCount));
	y_ = rcp (new mat_type (maxIterCount, 1));
      }
    // Even if the max iteration count hasn't changed, we still want
    // to fill in the projected least-squares problem data with zeros,
    // to ensure (for example) that R_ is upper triangular and that H_
    // is upper Hessenberg.
    (void) H_->putScalar(zero);
    (void) R_->putScalar(zero);
    (void) y_->putScalar(zero);

    lastUpdatedCol_ = -1; // column updates start with zero
    curNumIters_ = 0;
    flexible_ = needToAllocateZ;

    // This would be the place where subclasses may implement things
    // like orthogonalizing the first basis vector (after restart)
    // against the recycled subspace.
    postRestartHook();
  }

} // namespace Belos

#endif // __Belos_GmresBase_hpp
