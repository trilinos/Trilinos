// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file AnasaziTraceMin.hpp
  \brief Implementation of the trace minimization eigensolver
*/

#ifndef ANASAZI_TRACEMIN_HPP
#define ANASAZI_TRACEMIN_HPP

#include "AnasaziTypes.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziTraceMinBase.hpp"

#ifdef HAVE_ANASAZI_EPETRA
  #include "Epetra_Operator.h"
#endif

#include "AnasaziEigensolver.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "AnasaziMatOrthoManager.hpp"
#include "AnasaziSolverUtils.hpp"

#include "Teuchos_LAPACK.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseSolver.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"

// TODO: TraceMin performs some unnecessary computations upon restarting.  Fix it!

namespace Anasazi {
namespace Experimental {
/*! \class TraceMin

    \brief This class implements a TraceMIN iteration, a preconditioned iteration
    for solving linear symmetric positive definite eigenproblems.

    %TraceMin works by solving a constrained minimization problem
    \f[ \textrm{min trace}\left(Y^TKY\right) \textrm{ such that } Y^TMY = I\f]
    The Y which satisfies this problem is the set of eigenvectors corresponding to
    the eigenvalues of smallest magnitude.  While %TraceMin is capable of finding any
    eigenpairs of \f$KX = MX \Sigma\f$ (where K is symmetric and M symmetric positive
    definite), it converges to the eigenvalues of smallest magnitude of whatever it is
    given.  If a different set of eigenpairs is desired (such as the absolute smallest
    or the ones closest to a given shift), please perform a spectral transformation
    before passing the Eigenproblem to the solver.

    A %TraceMin iteration consists of the following operations:
    -# Solve the saddle point problem
          \f[\left[\begin{array}{lr} K & MX \\ X^TM & 0\end{array}\right]
          \left[\begin{array}{l} V \\ \tilde{L}\end{array}\right] =
          \left[\begin{array}{l} 0 \\ I\end{array}\right]\f]
    -# Form a section out of the current subspace V such that \f$X^TKX = \Sigma\f$ and
       \f$X^TMX = I\f$, where \f$\Sigma\f$ is diagonal.  \f$\left(\Sigma,X\right)\f$ is
       an approximation of the eigenpairs of smallest magnitude.
    -# Compute the new residual and check for convergence.

    The saddle point problem need not be solved to a low relative residual and can be
    solved either by directly forming the Schur complement or by using a projected
    Krylov subspace method to solve
    \f[ \left(PKP\right) \Delta = PMX \textrm{ with } P=I-BX\left(X^TB^2X\right)^{-1}X^TB\f]
    Then V is constructed as \f$V=X-\Delta\f$.  If a preconditioner H is used with the
    projected Krylov method, it is applied as
    \f[F=\left[I-H^{-1}BX\left(X^TBH^{-1}BX\right)^{-1}X^TB\right]H^{-1}\f]
    and we solve \f$FK\Delta = FMX\f$.
    (See A Parallel Implementation of the Trace Minimization Eigensolver by Eloy Romero
    and Jose E. Roman.)

    The convergence rate of %TraceMin is based on the distribution of eigenvalues.  If
    the eigenvalues are clustered far away from the origin, we have a slow rate of
    convergence.  We can improve our convergence rate by taking advantage of Ritz
    shifts.  Instead of solving \f$Kx=\lambda Mx\f$, we solve
    \f$\left(K-\sigma M\right) x = \left(\lambda - \sigma\right)Mx\f$.

    This method is described in <em>A Trace Minimization Algorithm for
    the Generalized Eigenvalue Problem</em>, Ahmed H. Sameh, John A.
    Wisniewski, SIAM Journal on Numerical Analysis, 19(6), pp. 1243-1259
    (1982)

    \ingroup anasazi_solver_framework

    \author Alicia Klinvex
*/
  template <class ScalarType, class MV, class OP>
  class TraceMin : public TraceMinBase<ScalarType,MV,OP> {
  public:
    //! @name Constructor/Destructor
    //@{

    /*! \brief %TraceMin constructor with eigenproblem, solver utilities, and
     *  parameter list of solver options.
     *
     * This constructor takes pointers required by the eigensolver, in addition
     * to a parameter list of options for the eigensolver. These options include
     * the following:
     *   - \c "Block Size" - an \c int specifying the subspace dimension used
     *     by the algorithm. This can also be specified using the setBlockSize()
     *     method.
     *   - \c "Maximum Iterations" - an \c int specifying the maximum number of
     *     %TraceMin iterations.
     *   - \c "Saddle Solver Type" - a \c string specifying the algorithm to use
     *     in solving the saddle point problem: "Schur Complement" or "Projected
     *     Krylov". Default: "Projected Krylov"
     *      - \c "Schur Complement": We explicitly form the Schur complement and
     *        use it to solve the linear system.  This option may be faster, but
     *        it is less numerically stable and does not ensure orthogonality
     *        between the current iterate and the update.
     *      - \c "Projected Krylov": Use a projected iterative method to solve
     *        the linear system. If %TraceMin was not given a preconditioner, it
     *        will use MINRES.  Otherwise, it will use GMRES.
     *   - \c "Shift Type" - a \c string specifying how to choose Ritz shifts:
     *        "No Shift", "Locked Shift", "Trace Leveled Shift", or "Original Shift".
     *        Default: "Trace Leveled Shift"
     *      - \c "No Shift": Do not perform Ritz shifts.  This option produces
     *        guaranteed convergence but converges linearly.  Not recommended.
     *      - \c "Locked Shift": Do not perform Ritz shifts until an eigenpair is
     *        locked.  Then, shift based on the largest converged eigenvalue.
     *        This option is roughly as safe as "No Shift" but may be faster.
     *      - \c "Trace Leveled Shift": Do not perform Ritz shifts until the trace
     *        of \f$X^TKX\f$ (i.e. the quantity being minimized has stagnated.
     *        Then, shift based on the strategy proposed in <em>The trace
     *        minimization method for the symmetric generalized eigenvalue problem.</em>
     *        Highly recommended.
     *      - \c "Original Shift": The original shifting strategy proposed in
     *        "The trace minimization method for the symmetric generalized
     *        eigenvalue problem."  Compute shifts based on the Ritz values,
     *        residuals, and clustering of the eigenvalues.  May produce incorrect
     *        results for indefinite matrices or small subspace dimensions.
     */
    TraceMin( const Teuchos::RCP<Eigenproblem<ScalarType,MV,OP> >    &problem,
                   const Teuchos::RCP<SortManager<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> > &sorter,
                   const Teuchos::RCP<OutputManager<ScalarType> >         &printer,
                   const Teuchos::RCP<StatusTest<ScalarType,MV,OP> >      &tester,
                   const Teuchos::RCP<MatOrthoManager<ScalarType,MV,OP> > &ortho,
                   Teuchos::ParameterList &params
                 );

  private:
    //
    // Convenience typedefs
    //
    typedef SolverUtils<ScalarType,MV,OP> Utils;
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename SCT::magnitudeType MagnitudeType;
    const MagnitudeType ONE;
    const MagnitudeType ZERO;
    const MagnitudeType NANVAL;

    // TraceMin specific methods
    void addToBasis(const Teuchos::RCP<const MV> Delta);

    void harmonicAddToBasis(const Teuchos::RCP<const MV> Delta);
  };

  //////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Implementations
  //
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor
  template <class ScalarType, class MV, class OP>
  TraceMin<ScalarType,MV,OP>::TraceMin(
        const Teuchos::RCP<Eigenproblem<ScalarType,MV,OP> >    &problem,
        const Teuchos::RCP<SortManager<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> > &sorter,
        const Teuchos::RCP<OutputManager<ScalarType> >         &printer,
        const Teuchos::RCP<StatusTest<ScalarType,MV,OP> >      &tester,
        const Teuchos::RCP<MatOrthoManager<ScalarType,MV,OP> > &ortho,
        Teuchos::ParameterList &params
        ) :
    TraceMinBase<ScalarType,MV,OP>(problem,sorter,printer,tester,ortho,params),
    ONE(Teuchos::ScalarTraits<MagnitudeType>::one()),
    ZERO(Teuchos::ScalarTraits<MagnitudeType>::zero()),
    NANVAL(Teuchos::ScalarTraits<MagnitudeType>::nan())
  {
  }


  template <class ScalarType, class MV, class OP>
  void TraceMin<ScalarType,MV,OP>::addToBasis(const Teuchos::RCP<const MV> Delta)
  {
    MVT::MvAddMv(ONE,*this->X_,-ONE,*Delta,*this->V_);

    if(this->hasM_)
    {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor lcltimer( *this->timerMOp_ );
#endif
      this->count_ApplyM_+= this->blockSize_;

      OPT::Apply(*this->MOp_,*this->V_,*this->MV_);
    }

    int rank;
    {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor lcltimer( *this->timerOrtho_ );
#endif

      if(this->numAuxVecs_ > 0)
      {
        rank = this->orthman_->projectAndNormalizeMat(*this->V_,this->auxVecs_,
               Teuchos::tuple(Teuchos::RCP< Teuchos::SerialDenseMatrix< int, ScalarType > >(Teuchos::null)),
               Teuchos::null,this->MV_,this->MauxVecs_);
      }
      else
      {
        rank = this->orthman_->normalizeMat(*this->V_,Teuchos::null,this->MV_);
      }
    }

    // FIXME (mfh 07 Oct 2014) This variable is currently unused, but
    // it would make sense to use it to check whether the block is
    // rank deficient.
    (void) rank;

    if(this->Op_ != Teuchos::null)
    {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor lcltimer( *this->timerOp_ );
#endif
      this->count_ApplyOp_+= this->blockSize_;
      OPT::Apply(*this->Op_,*this->V_,*this->KV_);
    }
  }



  template <class ScalarType, class MV, class OP>
  void TraceMin<ScalarType,MV,OP>::harmonicAddToBasis(const Teuchos::RCP<const MV> Delta)
  {
    // V = X - Delta
    MVT::MvAddMv(ONE,*this->X_,-ONE,*Delta,*this->V_);

    // Project out auxVecs
    if(this->numAuxVecs_ > 0)
    {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor lcltimer( *this->timerOrtho_ );
#endif
      this->orthman_->projectMat(*this->V_,this->auxVecs_);
    }

    // Compute KV
    if(this->Op_ != Teuchos::null)
    {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor lcltimer( *this->timerOp_ );
#endif
      this->count_ApplyOp_+= this->blockSize_;

      OPT::Apply(*this->Op_,*this->V_,*this->KV_);
    }

    // Normalize lclKV
    RCP< Teuchos::SerialDenseMatrix<int,ScalarType> > gamma = rcp(new Teuchos::SerialDenseMatrix<int,ScalarType>(this->blockSize_,this->blockSize_));
    int rank = this->orthman_->normalizeMat(*this->KV_,gamma);
    // FIXME (mfh 18 Feb 2015) It would make sense to check the rank.
    (void) rank;

    // lclV = lclV/gamma
    Teuchos::SerialDenseSolver<int,ScalarType> SDsolver;
    SDsolver.setMatrix(gamma);
    SDsolver.invert();
    RCP<MV> tempMV = MVT::CloneCopy(*this->V_);
    MVT::MvTimesMatAddMv(ONE,*tempMV,*gamma,ZERO,*this->V_);

    // Compute MV
    if(this->hasM_)
    {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor lcltimer( *this->timerMOp_ );
#endif
      this->count_ApplyM_+= this->blockSize_;

      OPT::Apply(*this->MOp_,*this->V_,*this->MV_);
    }
  }

}} // End of namespace Anasazi

#endif

// End of file AnasaziTraceMin.hpp
