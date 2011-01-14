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
#include <Teuchos_SerialDenseVector.hpp>
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
    GmresCantExtendBasis(const std::string& what_arg) : BelosError(what_arg) {}
  };

  /// \brief Candidate "basis" isn't a basis
  ///
  /// Thrown by GmresBase::advance(), when it rejects the computed
  /// candidate basis vector(s) due to (numerical) rank deficiency, 
  /// and doesn't know how to recover.
  ///
  /// This usually means that the candidate basis vectors from
  /// extendBasis() are not full rank after orthogonalization.
  /// CA-GMRES may choose to retry with a shorter candidate basis
  /// length, but if the candidate basis length is too short, it may
  /// opt to "give up."  In that case, advance() throws this
  /// exception.  Restarting with standard GMRES may be a good idea in
  /// that case.
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
    GmresRejectsCandidateBasis(const std::string& what_arg) : BelosError(what_arg) {}
  };


  /// \class GmresBase
  /// \brief Common state and functionality for new GMRES implementations
  /// \author Mark Hoemmen
  ///
  /// This class includes both state and functionality that are useful
  /// for different implementations of GMRES.  It does not implement
  /// the actual iterations; this is left for subclasses.
  /// Furthermore, it does not implement features like recycling, but
  /// it does include hooks for subclasses to add in such features.
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
    GmresBase (const Teuchos::RCP< LinearProblem<Scalar, MV, OP> >& problem,
	       const Teuchos::RCP< const OrthoManager<Scalar, MV> >& ortho,
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
    magnitude_type currentNativeResidualNorm ();

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
    magnitude_type arnoldiRelationError () const;

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
    virtual bool canExtendBasis() = 0;

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
    virtual void 
    updateUpperHessenbergMatrix (const Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> >& C_V,
				 const Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> >& B_V,
				 const Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> >& C_Z,
				 const Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> >& B_Z) = 0;

    //! Whether the subclass "accepted" the candidate basis
    virtual bool acceptedCandidateBasis() = 0;

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
			   const Teuchos::RCP<mat_type>& C_Z
			   const Teuchos::RCP<mat_type>& B_Z) {}
    //@}

  private:

    //! \name Helper methods
    //@{

    /// \brief Set first basis vector and initial residual norm
    ///
    /// \param r0 [in] Initial residual vector (left preconditioned,
    ///   if applicable)
    ///
    /// \note This is a helper method, to be called by the constructor
    ///   and by restart() only.
    void setInitialResidual(const Teuchos::RCP<const MV>& r0);

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
    magnitude_type updateProjectedLeastSquaresProblem ();

    /// Views of previously orthogonalized basis vectors
    ///
    /// \param V_prv [out] View of previously orthogonalized V basis
    ///   vectors (a view into V_)
    ///
    /// \param Z_prv [out] If running Flexible GMRES, view of
    ///   previously orthogonalized Z basis vectors (a view into Z_).
    ///   Else, Teuchos::null.
    ///
    /// \warning Don't call this until after setInitialResidual() has
    ///   been called (it's called in the constructor, and in
    ///   restart()).
    void 
    previousVectors (Teuchos::RCP<const MV>& V_prv,
		     Teuchos::RCP<const MV>& Z_prv) const;

    /// \brief Initial residual vector
    ///
    /// Initial residual vector (left preconditioned, if there is a
    /// left preconditioner) from the linear problem to solve.  This
    /// changes if the linear problem is updated.
    ///
    /// \note We keep this here so that it's easy to generate new
    /// basis vectors, since the V_ basis comes from the same space as
    /// the initial residual vector.
    ///
    /// \warning Don't call this until after the linear problem (lp_)
    ///   has been set.
    Teuchos::RCP<const MV> initResVec() const;

    //@}

  protected:
    //! \name Member data, inherited by subclasses
    //@{

    //! The linear problem to solve
    Teuchos::RCP< LinearProblem<Scalar, MV, OP> > lp_;

    //! Orthogonalization manager
    Teuchos::RCP< const OrthoManager<Scalar, MV> > ortho_;

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
    /// Storage for the native residual vector is cached and
    /// reallocated only when the LinearProblem changes.  (This
    /// ensures that the dimensions and data distribution are correct,
    /// i.e., are the same as the initial residual vector in the
    /// linear problem).
    Teuchos::RCP<MV> nativeResVec_;
    
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
    Teuchos::RCP<Teuchos::SerialDenseVector<int, Scalar> > y_;

    /// \brief Current RHS of the projected least-squares problem
    ///
    /// The current right-hand side of the projected least-squares
    /// problem \fn$\min_y \|\underline{H} y - \beta e_1\|_2\fn$.  z_
    /// starts out as \fn$\beta e_1\fn$ (where \fn$\beta\fn$ is the
    /// initial residual norm).  It is updated progressively along
    /// with the QR factorization of H_.
    Teuchos::RCP<Teuchos::SerialDenseVector<int, Scalar> > z_;

    /// \brief The initial residual norm
    ///
    /// GMRES makes use of the initial residual norm for solving the
    /// projected least-squares problem for the solution update
    /// coefficients.  In that case, it is usually called
    /// \fn$\beta\fn$.  For left-preconditioned GMRES, this is the
    /// preconditioned initial residual norm, else it's the
    /// unpreconditioned version.
    magnitude_type initialResidualNorm_;

    /// Last column of H_ for which the QR factorization (implicitly
    /// stored in theCosines, theSines, and R_) has been computed.
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
    
    // getNumIters() does not include the initial basis vector.
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

    RCP<const MV> V_prv, Z_prv, V_cur, Z_cur;
    previousVectors (V_prv, Z_prv);
    extendBasis (V_cur, Z_cur);
    preOrthogonalizeHook (V_cur, Z_cur);
    const int s = MVT::GetNumberVecs (*V_cur);

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
    // vector(s), and perhaps also to facilitate calling advance()
    // recursively (e.g., setting the "last rejected rank" so that the
    // next call of advance() computes fewer candidate basis vectors).
    if (! acceptedCandidateBasis())
      {
	TEST_FOR_EXCEPTION( !shouldRetryAdvance(), 
			    GmresRejectsCandidateBasis, 
			    "GMRES rejects the computed candidate basis vector" 
			    << ((flexible_ || s > 1) ? "s" : "")
			    << "." );
	advance(); // advance() recursively
      }
    else
      {
	curNumIters_ += rank;
	// This is done lazily, whenever the "native" residual norm
	// or the current solution update are requested.
	//
	//updateProjectedLeastSquaresProblem();
      }
  }


  template<class Scalar, class MV, class OP>
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType
  GmresBase<Scalar, MV, OP>::arnoldiRelationError () const
  {
    using Teuchos::Range1D;
    using Teuchos::RCP;
    const Scalar one = Teuchos::ScalarTraits<Scalar>::one();
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
    magnitude_type result (0);
    for (k = 0; k < m; ++k)
      result += norms[k];
    return Teuchos::ScalarTraits<magnitude_type>::squareroot(result);
  }


  template<class Scalar, class MV, class OP>
  GmresBase<Scalar, MV, OP>::
  GmresBase (const Teuchos::RCP< LinearProblem<Scalar, MV, OP> >& problem,
	     const Teuchos::RCP< const OrthoManager<Scalar, MV> >& ortho,
	     const int maxIterCount,
	     const bool flexible)
    : lp_ (problem), ortho_ (ortho)
  {
    using Teuchos::is_null;
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::rcp;
    const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
    const bool haveLeftPrec = ! is_null (lp_->getLeftPrec());
    const bool haveRightPrec = ! is_null (lp_->getRightPrec());
    
    // The solution update is always in the same space as the initial
    // solution guess (x0).
    RCP<const MV> x0 = lp_->getLHS();
    TEST_FOR_EXCEPTION(MVT::GetNumberVecs (*x0) != 1, std::invalid_argument,
		       "Our Arnoldi/GMRES implementation only works for "
		       "single-vector problems, but the supplied initial "
		       "guess has " << MVT::GetNumberVecs(*x0) << " columns.");
    xUpdate_ = MVT::Clone (x0, 1);
    MVT::MvInit (*xUpdate_, zero);

    // (Left preconditioned, if applicable) "exact" residual vector.
    RCP<const MV> r0 = haveLeftPrec ? lp_->getInitPrecResVec() : lp_->getInitResVec();
    TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*r0) != 1, 
		       std::invalid_argument,
		       "Our Arnoldi/GMRES implementation only works for "
		       "single-vector problems, but the supplied initial "
		       "residual vector has " 
		       << MVT::GetNumberVecs(*r0_) << " columns.");
    // If left preconditioning is used, then the "native" residual
    // vector is in the same space as the preconditioned "exact"
    // residual vector.  Otherwise, it is in the same space as the
    // right-hand side b and the unpreconditioned "exact" residual
    // vector.
    nativeResVec_ = MVT::CloneCopy (r0);

    // If left preconditioning is used, then the V_ vectors are in the
    // same space as the preconditioned "exact" residual vector.
    // Otherwise they are in the same space as the right-hand side b
    // and the unpreconditioned "exact" residual vector.
    V_ = MVT::Clone (r0, maxIterCount+1);

    // The Z_ vectors, if we need them, are always in the same space
    // as the initial solution guess (x0).  Even if the caller asks
    // for the flexible option, we don't allocate space for Z_ unless
    // the caller specifies a right preconditioner.
    Z_ = (flexible && haveRightPrec) ? MVT::Clone(x0, maxIterCount+1) : null;

    // These matrices and vectors encode the small dense projected
    // least-squares problem.
    H_ = rcp (new mat_type (maxIterCount+1, maxIterCount));
    R_ = rcp (new mat_type (maxIterCount+1, maxIterCount));
    y_ = rcp (new vec_type (maxIterCount));
    z_ = rcp (new vec_type (maxIterCount+1));
    // These cosines and sines encode the Q factor in the QR
    // factorization of the upper Hessenberg matrix H_.  We compute
    // this by computing H_ into R_ and operating on R_ in place; H_
    // itself is left alone.
    theCosines_.size(maxIterCount);
    theSines_.size(maxIterCount);

    // Compute initial residual norm, and set the first column of V_
    // to the scaled initial residual vector.
    setInitialResidual(r0);

    lastUpdatedCol_ = -1; // column updates start with zero
    curNumIters_ = 0;
    maxNumIters_ = maxIterCount;
    flexible_ = (flexible && haveRightPrec);
  }

  template<class Scalar, class MV, class OP>
  Teuchos::RCP<MV> 
  GmresBase<Scalar, MV, OP>::currentNativeResidualVector ()
  {
    const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
    const Scalar one = Teuchos::ScalarTraits<Scalar>::one();
    using Teuchos::is_null;
    using Teuchos::Range1D;
    using Teuchos::RCP;

    if (is_null (nativeResVec_) || MVT::GetNumberVecs (*nativeResVec_) != MVT::GetNumberVecs (*initResVec()))
      nativeResVec_ = MVT::CloneCopy (*initResVec(), Range1D(0,0));
    else 
      // Assign initial residual vector to nativeResVec; we will
      // update this below if the number of iterations is > 0.
      //
      // TODO method doesn't exist for Range1D inputs yet
      MVT::SetBlock (*initResVec(), Range1D(0,0), nativeResVec_);

    if (curNumIters_ > 0)
      {
	TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*V_) < curNumIters_+1, std::logic_error,
			   "Only " << MVT::GetNumberVecs(*V_) << " basis vectors "
			   "were given, but curNumIters+1=" << (curNumIters_+1) 
			   << " of them are required.  "
			   "This likely indicates a bug in Belos.");
	TEST_FOR_EXCEPTION(H_->numRows() < curNumIters_+1, std::logic_error,
			   "H only has " << H_->numRows() << " rows, but " 
			   << (curNumIters_+1) << " rows are required.  "
			   "This likely indicates a bug in Belos.");
	TEST_FOR_EXCEPTION(H_->numCols() < curNumIters_, std::logic_error,
			   "H only has " << H_->numCols() << " columns, but "
			   << (curNumIters_+1) << " columns are required.  "
			   "This likely indicates a bug in Belos.");
	TEST_FOR_EXCEPTION(y_->length() < curNumIters_, std::logic_error,
			   "y only has " << y_->length() << " entries, but "
			   curNumIters_ << " entries are required.  "
			   "This likely indicates a bug in Belos.");
	RCP<const MV> V_view = MVT::CloneView (*V_, Range1D(0, curNumIters_));
	const mat_type H_view (Teuchos::View, *H_, curNumIters_+1, curNumIters_);
	// SerialDenseVector doesn't have a nice subview copy constructor.
	const vec_type y_view (Teuchos::View, y_->values(), curNumIters_);
	vec_type H_times_y (numIters+1, 1);
	{
	  const int err = 
	    H_times_y.multiply (Teuchos::NOTRANS, Teuchos::NOTRANS, 
				one, H_view, y_view, zero);
	  TEST_FOR_EXCEPTION(err != 0, std::logic_error,
			     "In GMRES, when computing the current native resi"
			     "dual vector via the Arnoldi relation, H*y failed"
			     " due to incompatible dimensions.  This is likely"
			     " a Belos bug.");
	}
	// nativeResVec_ := V * (H * y) - r0 (where we stored r0 in nativeResVec_)
	MVT::MvTimesMatAddMv (one, V, H_times_y, -one, nativeResVec_);
      }
    return nativeResVec_;
  }

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
    const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
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
      return Teuchos::ScalarTraits<Scalar>::magnitude (z[endCol+1]);

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
    vec_type& z = *z_;
    
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
	blas.ROT (1, &z[j], 1, &z[j+1], 1, &theCosine, &theSine);
      }
    // Remember the rightmost updated column.
    lastUpdatedCol_ = endCol;

    // The absolute value of the last element of z gives the current
    // "native" residual norm.
    const magnitude_type nativeResNorm = 
      Teuchos::ScalarTraits<scalar_type>::magnitude (z[endCol+1]);

    // Now that we have the updated R factor of H, and the updated
    // right-hand side z, solve the least-squares problem by solving
    // the linear system Ry=z.
    // 
    // TRSM() overwrites the right-hand side with the solution, so
    // copy z into y.
    {
      vec_type y_view (Teuchos::View, *y_, endCol+1, 1);
      vec_type z_view (Teuchos::View, *z_, endCol+1, 1);
      y_view.assign (z_view);
    }
    // Solve Ry = z for y.
    blas.TRSM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NOTRANS, 
	      Teuchos::NON_UNIT_DIAG, numIters, 1, 
	      one, R_->values(), LDR, y_->values(), y_->stride());
    return nativeResNorm;
  }


  template<class Scalar, class MV, class OP>
  std::pair< Teuchos::RCP<MV>, typename Teuchos::ScalarTraits<Scalar>::magnitudeType >
  GmresBase<Scalar, MV, OP>::getCurrentUpdate ()
  {
    using Teuchos::is_null;
    using Teuchos::Range1D;
    using Teuchos::RCP;
    const Scalar one = Teuchos::ScalarTraits<Scalar>::one();
    const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();

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
	  xUpdate_ = MVT::Clone (Z_, 1);
	else
	  xUpdate_ = MVT::Clone (V_, 1);
      }
    if (curNumIters_ == 0)
      {
	MVT::MvInit (*xUpdate_, zero);
	return;
      }

    // When one iteration of GMRES is completed, the upper
    // Hessenberg matrix H is 2 by 1.  Thus, we want to update
    // columns up to and including curNumIters_ - 1: actually,
    // columns [lastUpdatedCol_+1, curNumIters_-1] (inclusive).  The
    // solution update gives you the current "native" residual norm
    // for free, so we save and return that as well.
    const magnitude_type nativeResNorm = updateProjectedLeastSquaresProblem();

    // Current iteration count does not include the initial basis vector.
    const int m = getNumIters();
    // Basis for solution update coefficients; Flexible GMRES uses the
    // preconditioned basis here.
    RCP<const MV> Z_view = flexible_ ? 
      MVT::CloneView (*Z_, Range1D(0, m-1)) : 
      MVT::CloneView (*V_, Range1D(0, m-1));
    // SerialDenseVector doesn't have a proper copy view constructor.
    const vec_type y_view (Teuchos::View, y_->values(), m);
    // xUpdate_ := Z_view * y_view
    MVT::MvTimesMatAddMv (one, *Z_view, y_view, zero, *xUpdate_);
    return std::make_pair (xUpdate_, nativeResNorm);
  }

  template<class Scalar, class MV, class OP>
  void
  GmresBase<Scalar, MV, OP>::
  setInitialResidual(const Teuchos::RCP<const MV>& r0)
  {
    using Teuchos::Range1D;
    using Teuchos::RCP;

    // Initial residual norm is computed with respect to the inner
    // product defined by the OrthoManager.  Since the basis
    // vectors are orthonormal (resp. unitary) with respect to that
    // inner product, this ensures that the least-squares minimization
    // is also computed with respect to that inner product.
    mat_type result(1,1);
    ortho_->innerProd (*r0, *r0, result);
    // Norm is a magnitude_type.  innerProd() does not promise a
    // strictly nonzero imaginary part, but it should be "relatively
    // small."  Defining "relatively small" requires some error
    // analysis, which we omit; instead, we simply take the magnitude
    // (which hopefully is nonnegative by construction, if complex
    // absolute value was implemented correctly).
    initialResidualNorm_ = Teuchos::ScalarTraits<Scalar>::magnitude (result[0,0]);

    // Initialize right-hand side of projected least-squares problem
    (void) z_.putScalar(zero);
    z_[0] = Scalar (initialResidualNorm_);

    MVT::SetBlock (r0, Range1D(0,0), V_);
    RCP<MV> v1 = MVT::CloneViewNonConst (V_, Range1D(0,0));
    MVT::MvScale (v1, Scalar(1)/initialResidualNorm_);
  }

  template<class Scalar, class MV, class OP>
  void 
  GmresBase<Scalar, MV, OP>::backOut (const int numIters)
  {
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
    (*z_)[0] = Scalar (initialResidualNorm_);

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
    using Teuchos::is_null;
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::rcp;

    // This would be the place where subclasses may implement things
    // like recycling basis vectors for the next restart cycle.  Note
    // that the linear problem hasn't been updated yet; this happens
    // below.
    preRestartHook (lp_, maxIterCount);

    // Make sure that the LinearProblem object gets the current
    // approximate solution (which becomes the initial guess for the
    // upcoming new restart cycle), so that the "exact" residuals are
    // correct.
    (void) getCurrentUpdate (); // results stored in xUpdate_
    (void) lp_->updateSolution (xUpdate_, true, Scalar(1.0));

    const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
    const bool haveLeftPrec = ! is_null (lp_->getLeftPrec());
    const bool haveRightPrec = ! is_null (lp_->getRightPrec());
    // Initial residual vector (left preconditioned, if applicable).
    RCP<const MV> r0 = haveLeftPrec ? lp_->getInitPrecResVec() : lp_->getInitResVec();

    // If the new max iteration count has changed, reallocate space
    // for basis vectors and the projected least-squares problem.
    if (maxNumIters_ != maxIterCount)
      {
	maxNumIters_ = maxIterCount;
	V_ = MVT::Clone (r0, maxIterCount+1);
	Z_ = (flexible_ && haveRightPrec) ? MVT::Clone(lp.getLHS(), maxIterCount+1) : null;
	H_ = rcp (new mat_type (maxIterCount+1, maxIterCount));
	R_ = rcp (new mat_type (maxIterCount+1, maxIterCount));
	y_ = rcp (new vec_type (maxIterCount));
	z_ = rcp (new vec_type (maxIterCount+1));
      }
    // Even if the max iteration count hasn't changed, we still want
    // to fill in the projected least-squares problem data with zeros,
    // to ensure (for example) that R_ is upper triangular and that H_
    // is upper Hessenberg.
    (void) H_->putScalar(zero);
    (void) R_->putScalar(zero);
    (void) y_->putScalar(zero);
    (void) z_->putScalar(zero);

    setInitialResidual(r0);

    lastUpdatedCol_ = -1; // column updates start with zero
    curNumIters_ = 0;
    flexible_ = (flexible && haveRightPrec);

    // This would be the place where subclasses may implement things
    // like orthogonalizing the first basis vector (after restart)
    // against the recycled subspace.
    postRestartHook();
  }

  template<class Scalar, class MV, class OP>
  Teuchos::RCP<const MV> GmresBase<Scalar, MV, OP>::initResVec() const
  {
    const bool haveLeftPrec = ! is_null (lp_->getLeftPrec());
    const bool haveRightPrec = ! is_null (lp_->getRightPrec());
    // Initial residual vector (left preconditioned, if applicable).
    return = haveLeftPrec ? lp_->getInitPrecResVec() : lp_->getInitResVec();
  }

} // namespace Belos

#endif // __Belos_GmresBase_hpp
