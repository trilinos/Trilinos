// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//! \file AnasaziSIRTR.hpp


#ifndef ANASAZI_SIRTR_HPP
#define ANASAZI_SIRTR_HPP

#include "AnasaziTypes.hpp"
#include "AnasaziRTRBase.hpp"

#include "AnasaziEigensolver.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "Teuchos_LAPACK.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"

/*!     \class Anasazi::SIRTR

        SIRTR ("skinny IRTR") is a non-caching, lower-memory implementation of the Implicit
        Riemannian Trust-Region (%IRTR) eigensolver.

        The solver uses between 6 and 8 blocks of vectors, compared to the
        requirements by IRTR of 10 to 13 blocks of vectors. The base requirement
        is 6 blocks of vectors, where a block of vectors contains a number of vectors equal to the
        block size specified for the solver (see RTRBase::getBlockSize()).
        Additional blocks are required when solving a generalized eigenvalue problem or when using a preconditioiner.

        For more information, see RTRBase.

        \ingroup anasazi_solver_framework

        \author Chris Baker
*/


// TODO: add randomization
// TODO: add expensive debug checking on Teuchos_Debug

namespace Anasazi {

  template <class ScalarType, class MV, class OP>
  class SIRTR : public RTRBase<ScalarType,MV,OP> {
  public:

    //! @name Constructor/Destructor
    //@{

    /*! \brief %SIRTR constructor with eigenproblem, solver utilities, and parameter list of solver options.
     *
     * This constructor takes pointers required by the eigensolver, in addition
     * to a parameter list of options for the eigensolver. These options include the following:
     *   - "Rho Prime" - an \c MagnitudeType specifying the size of the implicit trust-region radius.
     *   - "Block Size" - an \c int specifying the block size used by the algorithm. This can also be specified using the setBlockSize() method.
     *   - "Leftmost" - a \c bool specifying whether the solver is computing the
     *     leftmost ("SR") or rightmost ("LR") eigenvalues. Default: true. This must be in accord with the SortManager pass to the constructor.
     *   - "Kappa Convergence" - a \c MagnitudeType specifing the rate of convergence for the linear convergence regime. Default: 0.1
     *   - "Theta Convergence" - a \c MagnitudeType specifing the order of convergence for the linear convergence regime. theta implies a convergence order of theta+1. Default: 1.0
     */
    SIRTR( const Teuchos::RCP<Eigenproblem<ScalarType,MV,OP> >    &problem,
           const Teuchos::RCP<SortManager<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> >           &sorter,
           const Teuchos::RCP<OutputManager<ScalarType> >         &printer,
           const Teuchos::RCP<StatusTest<ScalarType,MV,OP> >      &tester,
           const Teuchos::RCP<GenOrthoManager<ScalarType,MV,OP> > &ortho,
           Teuchos::ParameterList                                 &params
        );

    //! %SIRTR destructor
    virtual ~SIRTR() {};

    //@}

    //! @name Solver methods
    //@{

    //! \brief Impemements Eigensolver. The outer %IRTR iteration. See RTRBase::iterate().
    void iterate();

    //@}

    //!  @name Output methods
    //@{

    //! Impemements Eigensolver. This method requests that the solver print out its current status to screen.
    void currentStatus(std::ostream &os);

    //@}

  private:
    //
    // Convenience typedefs
    //
    typedef SolverUtils<ScalarType,MV,OP> Utils;
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename SCT::magnitudeType MagnitudeType;
    typedef Teuchos::ScalarTraits<MagnitudeType> MAT;
    enum trRetType {
      UNINITIALIZED = 0,
      MAXIMUM_ITERATIONS,
      NEGATIVE_CURVATURE,
      EXCEEDED_TR,
      KAPPA_CONVERGENCE,
      THETA_CONVERGENCE
    };
    // these correspond to above
    std::vector<std::string> stopReasons_;
    //
    // Consts
    //
    const MagnitudeType ZERO;
    const MagnitudeType ONE;
    //
    // Internal methods
    //
    //! \brief The inner %IRTR iteration. See RTRBase::solveTRSubproblem().
    void solveTRSubproblem();
    //
    // rho_prime
    MagnitudeType rho_prime_;
    //
    // norm of initial gradient: this is used for scaling
    MagnitudeType normgradf0_;
    //
    // tr stopping reason
    trRetType innerStop_;
    //
    // number of inner iterations
    int innerIters_, totalInnerIters_;
  };




  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor
  template <class ScalarType, class MV, class OP>
  SIRTR<ScalarType,MV,OP>::SIRTR(
        const Teuchos::RCP<Eigenproblem<ScalarType,MV,OP> >    &problem,
        const Teuchos::RCP<SortManager<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> >           &sorter,
        const Teuchos::RCP<OutputManager<ScalarType> >         &printer,
        const Teuchos::RCP<StatusTest<ScalarType,MV,OP> >      &tester,
        const Teuchos::RCP<GenOrthoManager<ScalarType,MV,OP> > &ortho,
        Teuchos::ParameterList                                 &params
        ) :
    RTRBase<ScalarType,MV,OP>(problem,sorter,printer,tester,ortho,params,"SIRTR",true),
    ZERO(MAT::zero()),
    ONE(MAT::one()),
    totalInnerIters_(0)
  {
    // set up array of stop reasons
    stopReasons_.push_back("n/a");
    stopReasons_.push_back("maximum iterations");
    stopReasons_.push_back("negative curvature");
    stopReasons_.push_back("exceeded TR");
    stopReasons_.push_back("kappa convergence");
    stopReasons_.push_back("theta convergence");

    rho_prime_ = params.get("Rho Prime",0.5);
    TEUCHOS_TEST_FOR_EXCEPTION(rho_prime_ <= 0 || rho_prime_ >= 1,std::invalid_argument,
                       "Anasazi::SIRTR::constructor: rho_prime must be in (0,1).");
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // TR subproblem solver
  //
  // FINISH:
  //   define pre- and post-conditions
  //
  // POST:
  //   delta_,Adelta_,Hdelta_ undefined
  //
  template <class ScalarType, class MV, class OP>
  void SIRTR<ScalarType,MV,OP>::solveTRSubproblem() {

    // return one of:
    // MAXIMUM_ITERATIONS
    // NEGATIVE_CURVATURE
    // EXCEEDED_TR
    // KAPPA_CONVERGENCE
    // THETA_CONVERGENCE

    using Teuchos::RCP;
    using Teuchos::tuple;
    using Teuchos::null;
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
    using Teuchos::TimeMonitor;
#endif
    using std::endl;
    typedef Teuchos::RCP<const MV> PCMV;
    typedef Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > PSDM;

    innerStop_ = MAXIMUM_ITERATIONS;

    const int n = MVT::GetGlobalLength(*this->eta_);
    const int p = this->blockSize_;
    const int d = n*p - (p*p+p)/2;

    // We have the following:
    //
    // X'*B*X = I
    // X'*A*X = theta_
    //
    // We desire to remain in the trust-region:
    // { eta : rho_Y(eta) \geq rho_prime }
    // where
    // rho_Y(eta) = 1/(1+eta'*B*eta)
    // Therefore, the trust-region is
    // { eta : eta'*B*eta \leq 1/rho_prime - 1 }
    //
    const double D2 = ONE/rho_prime_ - ONE;

    std::vector<MagnitudeType> d_Hd(p), alpha(p), beta(p), z_r(p), zold_rold(p);
    std::vector<MagnitudeType> eBe(p), eBd(p), dBd(p), new_eBe(p);
    MagnitudeType r0_norm;

    MVT::MvInit(*this->eta_ ,0.0);

    //
    // R_ contains direct residuals:
    //    R_ = A X_ - B X_ diag(theta_)
    //
    // r0 = grad f(X) = 2 P_BX A X = 2 P_BX (A X - B X diag(theta_) = 2 proj(R_)
    // We will do this in place.
    // For seeking the rightmost, we want to maximize f
    // This is the same as minimizing -f
    // Substitute all f with -f here. In particular,
    //    grad -f(X) = -grad f(X)
    //    Hess -f(X) = -Hess f(X)
    //
    {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
      TimeMonitor lcltimer( *this->timerOrtho_ );
#endif
      this->orthman_->projectGen(
          *this->R_,                                            // operating on R
          tuple<PCMV>(this->BV_),tuple<PCMV>(this->V_),false,   // P_{BV,V}, and <BV,V>_B != I
          tuple<PSDM>(null),                                    // don't care about coeffs
          null,tuple<PCMV>(null), tuple<PCMV>(this->BV_));      // don't have B*BV, but do have B*V
      if (this->leftMost_) {
        MVT::MvScale(*this->R_,2.0);
      }
      else {
        MVT::MvScale(*this->R_,-2.0);
      }
    }

    r0_norm = MAT::squareroot( RTRBase<ScalarType,MV,OP>::ginner(*this->R_) );
    //
    // kappa (linear) convergence
    // theta (superlinear) convergence
    //
    MagnitudeType kconv = r0_norm * this->conv_kappa_;
    // FINISH: consider inserting some scaling here
    // MagnitudeType tconv = r0_norm * MAT::pow(r0_norm/normgradf0_,this->conv_theta_);
    MagnitudeType tconv = MAT::pow(r0_norm,this->conv_theta_+ONE);
    if (this->om_->isVerbosity(Debug)) {
      this->om_->stream(Debug)
        << " >> |r0|       : " << r0_norm << endl
        << " >> kappa conv : " << kconv << endl
        << " >> theta conv : " << tconv << endl;
    }

    //
    // For Olsen preconditioning, the preconditioner is
    // Z = P_{Prec^-1 BX, BX} Prec^-1 R
    // for efficiency, we compute Prec^-1 BX once here for use later
    // Otherwise, we don't need PBX
    if (this->hasPrec_ && this->olsenPrec_)
    {
      std::vector<int> ind(this->blockSize_);
      for (int i=0; i<this->blockSize_; ++i) ind[i] = this->numAuxVecs_+i;
      Teuchos::RCP<MV> PBX = MVT::CloneViewNonConst(*this->PBV_,ind);
      {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        TimeMonitor prectimer( *this->timerPrec_ );
#endif
        OPT::Apply(*this->Prec_,*this->BX_,*PBX);
        this->counterPrec_ += this->blockSize_;
      }
      PBX = Teuchos::null;
    }

    // Z = P_{Prec^-1 BV, BV} Prec^-1 R
    //    Prec^-1 BV in PBV
    // or
    // Z = P_{BV,BV} Prec^-1 R
    if (this->hasPrec_)
    {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
      TimeMonitor prectimer( *this->timerPrec_ );
#endif
      OPT::Apply(*this->Prec_,*this->R_,*this->Z_);
      this->counterPrec_ += this->blockSize_;
      // the orthogonalization time counts under Ortho and under Preconditioning
      if (this->olsenPrec_) {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        TimeMonitor orthtimer( *this->timerOrtho_ );
#endif
        this->orthman_->projectGen(
            *this->Z_,                                             // operating on Z
            tuple<PCMV>(this->PBV_),tuple<PCMV>(this->V_),false,   // P_{PBV,V}, B inner product, and <PBV,V>_B != I
            tuple<PSDM>(null),                                     // don't care about coeffs
            null,tuple<PCMV>(null), tuple<PCMV>(this->BV_));       // don't have B*PBV or B*Z, but do have B*V
      }
      else {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        TimeMonitor orthtimer( *this->timerOrtho_ );
#endif
        this->orthman_->projectGen(
            *this->Z_,                                             // operating on Z
            tuple<PCMV>(this->BV_),tuple<PCMV>(this->V_),false,    // P_{BV,V}, and <BV,V>_B != I
            tuple<PSDM>(null),                                     // don't care about coeffs
            null,tuple<PCMV>(null), tuple<PCMV>(this->BV_));       // don't have B*BV, but do have B*V
      }
      RTRBase<ScalarType,MV,OP>::ginnersep(*this->R_,*this->Z_,z_r);
    }
    else {
      // Z = R
      MVT::MvAddMv(ONE,*this->R_,ZERO,*this->R_,*this->Z_);
      RTRBase<ScalarType,MV,OP>::ginnersep(*this->R_,z_r);
    }

    if (this->om_->isVerbosity( Debug )) {
      // Check that gradient is B-orthogonal to X
      typename RTRBase<ScalarType,MV,OP>::CheckList chk;
      chk.checkBR = true;
      if (this->hasPrec_) chk.checkZ  = true;
      this->om_->print( Debug, this->accuracyCheck(chk, "after computing gradient") );
    }
    else if (this->om_->isVerbosity( OrthoDetails )) {
      // Check that gradient is B-orthogonal to X
      typename RTRBase<ScalarType,MV,OP>::CheckList chk;
      chk.checkBR = true;
      if (this->hasPrec_) chk.checkZ  = true;
      this->om_->print( OrthoDetails, this->accuracyCheck(chk, "after computing gradient") );
    }

    // delta = -z
    MVT::MvAddMv(-ONE,*this->Z_,ZERO,*this->Z_,*this->delta_);

    if (this->om_->isVerbosity(Debug)) {
      // compute the model at eta
      // we need Heta, which requires A*eta and B*eta
      // we also need A*X
      // use Z for storage of these
      std::vector<MagnitudeType> eAx(this->blockSize_),
        d_eAe(this->blockSize_),
        d_eBe(this->blockSize_),
        d_mxe(this->blockSize_);
      // compute AX and <eta,AX>
      {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        TimeMonitor lcltimer( *this->timerAOp_ );
#endif
        OPT::Apply(*this->AOp_,*this->X_,*this->Z_);
        this->counterAOp_ += this->blockSize_;
      }
      RTRBase<ScalarType,MV,OP>::ginnersep(*this->eta_,*this->Z_,eAx);
      // compute A*eta and <eta,A*eta>
      {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        TimeMonitor lcltimer( *this->timerAOp_ );
#endif
        OPT::Apply(*this->AOp_,*this->eta_,*this->Z_);
        this->counterAOp_ += this->blockSize_;
      }
      RTRBase<ScalarType,MV,OP>::ginnersep(*this->eta_,*this->Z_,d_eAe);
      // compute B*eta and <eta,B*eta>
      if (this->hasBOp_) {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        TimeMonitor lcltimer( *this->timerBOp_ );
#endif
        OPT::Apply(*this->BOp_,*this->eta_,*this->Z_);
        this->counterBOp_ += this->blockSize_;
      }
      else {
        MVT::MvAddMv(ONE,*this->eta_,ZERO,*this->eta_,*this->Z_);
      }
      RTRBase<ScalarType,MV,OP>::ginnersep(*this->eta_,*this->Z_,d_eBe);
      // compute model:
      //    m_x(eta) = theta + 2*eta'*A*x + eta'*A*eta - eta'*B*eta*theta
      if (this->leftMost_) {
        for (int j=0; j<this->blockSize_; ++j) {
          d_mxe[j] = this->theta_[j] + 2*eAx[j] + d_eAe[j] - d_eBe[j]*this->theta_[j];
        }
      }
      else {
        for (int j=0; j<this->blockSize_; ++j) {
          d_mxe[j] = -this->theta_[j] - 2*eAx[j] - d_eAe[j] + d_eBe[j]*this->theta_[j];
        }
      }
      this->om_->stream(Debug)
        << " Debugging checks: SIRTR inner iteration " << innerIters_ << endl
        << " >> m_X(eta) : " << std::accumulate(d_mxe.begin(),d_mxe.end(),0.0) << endl;
      for (int j=0; j<this->blockSize_; ++j) {
        this->om_->stream(Debug)
          << " >> m_X(eta_" << j << ") : " << d_mxe[j] << endl;
      }
    }

    ////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////
    // the inner/tCG loop
    for (innerIters_=1; innerIters_<=d; ++innerIters_) {

      //
      // [Hdelta,Adelta,Bdelta] = Hess*delta = 2 Proj(A*delta - B*delta*X'*A*X)
      // X'*A*X = diag(theta), so that
      // (B*delta)*diag(theta) can be done on the cheap
      //
      {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        TimeMonitor lcltimer( *this->timerAOp_ );
#endif
        OPT::Apply(*this->AOp_,*this->delta_,*this->Z_);
        this->counterAOp_ += this->blockSize_;
      }
      if (this->hasBOp_) {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        TimeMonitor lcltimer( *this->timerBOp_ );
#endif
        OPT::Apply(*this->BOp_,*this->delta_,*this->Hdelta_);
        this->counterBOp_ += this->blockSize_;
      }
      else {
        MVT::MvAddMv(ONE,*this->delta_,ZERO,*this->delta_,*this->Hdelta_);
      }
      // while we have B*delta, compute <eta,B*delta> and <delta,B*delta>
      // these will be needed below
      RTRBase<ScalarType,MV,OP>::ginnersep(*this->eta_  ,*this->Hdelta_,eBd);
      RTRBase<ScalarType,MV,OP>::ginnersep(*this->delta_,*this->Hdelta_,dBd);
      // put 2*A*d - 2*B*d*theta --> Hd
      {
        std::vector<ScalarType> theta_comp(this->theta_.begin(),this->theta_.end());
        MVT::MvScale(*this->Hdelta_,theta_comp);
      }
      if (this->leftMost_) {
        MVT::MvAddMv( 2.0,*this->Z_,-2.0,*this->Hdelta_,*this->Hdelta_);
      }
      else {
        MVT::MvAddMv(-2.0,*this->Z_, 2.0,*this->Hdelta_,*this->Hdelta_);
      }
      // apply projector
      {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        TimeMonitor lcltimer( *this->timerOrtho_ );
#endif
        this->orthman_->projectGen(
            *this->Hdelta_,                                       // operating on Hdelta
            tuple<PCMV>(this->BV_),tuple<PCMV>(this->V_),false,   // P_{BV,V}, and <BV,V>_B != I
            tuple<PSDM>(null),                                    // don't care about coeffs
            null,tuple<PCMV>(null), tuple<PCMV>(this->BV_));      // don't have B*BV, but do have B*V
      }
      RTRBase<ScalarType,MV,OP>::ginnersep(*this->delta_,*this->Hdelta_,d_Hd);


      // compute update step
      for (unsigned int j=0; j<alpha.size(); ++j)
      {
        alpha[j] = z_r[j]/d_Hd[j];
      }

      // compute new B-norms
      for (unsigned int j=0; j<alpha.size(); ++j)
      {
        new_eBe[j] = eBe[j] + 2*alpha[j]*eBd[j] + alpha[j]*alpha[j]*dBd[j];
      }

      if (this->om_->isVerbosity(Debug)) {
        for (unsigned int j=0; j<alpha.size(); j++) {
          this->om_->stream(Debug)
            << "     >> z_r[" << j << "]  : " << z_r[j]
            << "    d_Hd[" << j << "]  : " << d_Hd[j] << endl
            << "     >> eBe[" << j << "]  : " << eBe[j]
            << "    neweBe[" << j << "]  : " << new_eBe[j] << endl
            << "     >> eBd[" << j << "]  : " << eBd[j]
            << "    dBd[" << j << "]  : " << dBd[j] << endl;
        }
      }

      // check truncation criteria: negative curvature or exceeded trust-region
      std::vector<int> trncstep;
      trncstep.reserve(p);
      // trncstep will contain truncated step, due to
      //   negative curvature or exceeding implicit trust-region
      bool atleastonenegcur = false;
      for (unsigned int j=0; j<d_Hd.size(); ++j) {
        if (d_Hd[j] <= 0) {
          trncstep.push_back(j);
          atleastonenegcur = true;
        }
        else if (new_eBe[j] > D2) {
          trncstep.push_back(j);
        }
      }

      if (!trncstep.empty())
      {
        // compute step to edge of trust-region, for trncstep vectors
        if (this->om_->isVerbosity(Debug)) {
          for (unsigned int j=0; j<trncstep.size(); ++j) {
            this->om_->stream(Debug)
              << " >> alpha[" << trncstep[j] << "]  : " << alpha[trncstep[j]] << endl;
          }
        }
        for (unsigned int j=0; j<trncstep.size(); ++j) {
          int jj = trncstep[j];
          alpha[jj] = ( -eBd[jj] + MAT::squareroot(eBd[jj]*eBd[jj] + dBd[jj]*(D2-eBe[jj]) ) ) / dBd[jj];
        }
        if (this->om_->isVerbosity(Debug)) {
          for (unsigned int j=0; j<trncstep.size(); ++j) {
            this->om_->stream(Debug)
              << " >> tau[" << trncstep[j] << "]  : " << alpha[trncstep[j]] << endl;
          }
        }
        if (atleastonenegcur) {
          innerStop_ = NEGATIVE_CURVATURE;
        }
        else {
          innerStop_ = EXCEEDED_TR;
        }
      }

      // compute new eta = eta + alpha*delta
      // we need delta*diag(alpha)
      // do this in situ in delta_ and friends (we will note this for below)
      // then set eta_ = eta_ + delta_
      {
        std::vector<ScalarType> alpha_comp(alpha.begin(),alpha.end());
        MVT::MvScale(*this->delta_,alpha_comp);
        MVT::MvScale(*this->Hdelta_,alpha_comp);
      }
      MVT::MvAddMv(ONE,*this->delta_ ,ONE,*this->eta_ ,*this->eta_);

      // store new eBe
      eBe = new_eBe;

      //
      // print some debugging info
      if (this->om_->isVerbosity(Debug)) {
        // compute the model at eta
        // we need Heta, which requires A*eta and B*eta
        // we also need A*X
        // use Z for storage of these
        std::vector<MagnitudeType> eAx(this->blockSize_),
          d_eAe(this->blockSize_),
          d_eBe(this->blockSize_),
          d_mxe(this->blockSize_);
        // compute AX and <eta,AX>
        {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
          TimeMonitor lcltimer( *this->timerAOp_ );
#endif
          OPT::Apply(*this->AOp_,*this->X_,*this->Z_);
          this->counterAOp_ += this->blockSize_;
        }
        RTRBase<ScalarType,MV,OP>::ginnersep(*this->eta_,*this->Z_,eAx);
        // compute A*eta and <eta,A*eta>
        {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
          TimeMonitor lcltimer( *this->timerAOp_ );
#endif
          OPT::Apply(*this->AOp_,*this->eta_,*this->Z_);
          this->counterAOp_ += this->blockSize_;
        }
        RTRBase<ScalarType,MV,OP>::ginnersep(*this->eta_,*this->Z_,d_eAe);
        // compute B*eta and <eta,B*eta>
        if (this->hasBOp_) {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
          TimeMonitor lcltimer( *this->timerBOp_ );
#endif
          OPT::Apply(*this->BOp_,*this->eta_,*this->Z_);
          this->counterBOp_ += this->blockSize_;
        }
        else {
          MVT::MvAddMv(ONE,*this->eta_,ZERO,*this->eta_,*this->Z_);
        }
        RTRBase<ScalarType,MV,OP>::ginnersep(*this->eta_,*this->Z_,d_eBe);
        // compute model:
        //    m_x(eta) = theta + 2*eta'*A*x + eta'*A*eta - eta'*B*eta*theta
        if (this->leftMost_) {
          for (int j=0; j<this->blockSize_; ++j) {
            d_mxe[j] = this->theta_[j] + 2*eAx[j] + d_eAe[j] - d_eBe[j]*this->theta_[j];
          }
        }
        else {
          for (int j=0; j<this->blockSize_; ++j) {
            d_mxe[j] = -this->theta_[j] - 2*eAx[j] - d_eAe[j] + d_eBe[j]*this->theta_[j];
          }
        }
        this->om_->stream(Debug)
          << " Debugging checks: SIRTR inner iteration " << innerIters_ << endl
          << " >> m_X(eta) : " << std::accumulate(d_mxe.begin(),d_mxe.end(),0.0) << endl;
        for (int j=0; j<this->blockSize_; ++j) {
          this->om_->stream(Debug)
            << " >> m_X(eta_" << j << ") : " << d_mxe[j] << endl;
        }
      }

      //
      // if we found negative curvature or exceeded trust-region, then quit
      if (!trncstep.empty()) {
        break;
      }

      // update gradient of m
      // R = R + Hdelta*diag(alpha)
      // however, Hdelta_ already stores Hdelta*diag(alpha)
      // so just add them
      MVT::MvAddMv(ONE,*this->Hdelta_,ONE,*this->R_,*this->R_);
      {
        // re-tangentialize r
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        TimeMonitor lcltimer( *this->timerOrtho_ );
#endif
        this->orthman_->projectGen(
            *this->R_,                                            // operating on R
            tuple<PCMV>(this->BV_),tuple<PCMV>(this->V_),false,   // P_{BV,V}, and <BV,V>_B != I
            tuple<PSDM>(null),                                    // don't care about coeffs
            null,tuple<PCMV>(null), tuple<PCMV>(this->BV_));      // don't have B*BV, but do have B*V
      }

      //
      // check convergence
      MagnitudeType r_norm = MAT::squareroot(RTRBase<ScalarType,MV,OP>::ginner(*this->R_,*this->R_));

      //
      // check local convergece
      //
      // kappa (linear) convergence
      // theta (superlinear) convergence
      //
      if (this->om_->isVerbosity(Debug)) {
        this->om_->stream(Debug)
          << " >> |r" << innerIters_ << "|        : " << r_norm << endl;
      }
      if ( r_norm <= ANASAZI_MIN(tconv,kconv) ) {
        if (tconv <= kconv) {
          innerStop_ = THETA_CONVERGENCE;
        }
        else {
          innerStop_ = KAPPA_CONVERGENCE;
        }
        break;
      }

      // Z = P_{Prec^-1 BV, BV} Prec^-1 R
      //    Prec^-1 BV in PBV
      // or
      // Z = P_{BV,BV} Prec^-1 R
      zold_rold = z_r;
      if (this->hasPrec_)
      {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        TimeMonitor prectimer( *this->timerPrec_ );
#endif
        OPT::Apply(*this->Prec_,*this->R_,*this->Z_);
        this->counterPrec_ += this->blockSize_;
        // the orthogonalization time counts under Ortho and under Preconditioning
        if (this->olsenPrec_) {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
          TimeMonitor orthtimer( *this->timerOrtho_ );
#endif
          this->orthman_->projectGen(
              *this->Z_,                                             // operating on Z
              tuple<PCMV>(this->PBV_),tuple<PCMV>(this->V_),false,   // P_{PBV,V}, B inner product, and <PBV,V>_B != I
              tuple<PSDM>(null),                                     // don't care about coeffs
              null,tuple<PCMV>(null), tuple<PCMV>(this->BV_));       // don't have B*PBV or B*Z, but do have B*V
        }
        else {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
          TimeMonitor orthtimer( *this->timerOrtho_ );
#endif
          this->orthman_->projectGen(
              *this->Z_,                                             // operating on Z
              tuple<PCMV>(this->BV_),tuple<PCMV>(this->V_),false,    // P_{BV,V}, and <BV,V>_B != I
              tuple<PSDM>(null),                                     // don't care about coeffs
              null,tuple<PCMV>(null), tuple<PCMV>(this->BV_));       // don't have B*BV, but do have B*V
        }
        RTRBase<ScalarType,MV,OP>::ginnersep(*this->R_,*this->Z_,z_r);
      }
      else {
        // Z = R
        MVT::MvAddMv(ONE,*this->R_,ZERO,*this->R_,*this->Z_);
        RTRBase<ScalarType,MV,OP>::ginnersep(*this->R_,z_r);
      }

      // compute new search direction
      // below, we need to perform
      //   delta = -Z + delta*diag(beta)
      // however, delta_ currently stores delta*diag(alpha)
      // therefore, set
      //   beta_ to beta/alpha
      // so that
      //   delta_ = delta_*diag(beta_)
      // will in fact result in
      //   delta_ = delta_*diag(beta_)
      //          = delta*diag(alpha)*diag(beta/alpha)
      //          = delta*diag(beta)
      // i hope this is numerically sound...
      for (unsigned int j=0; j<beta.size(); ++j) {
        beta[j] = z_r[j]/(zold_rold[j]*alpha[j]);
      }
      {
        std::vector<ScalarType> beta_comp(beta.begin(),beta.end());
        MVT::MvScale(*this->delta_,beta_comp);
      }
      MVT::MvAddMv(-ONE,*this->Z_,ONE,*this->delta_,*this->delta_);

    }
    // end of the inner iteration loop
    //////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////
    if (innerIters_ > d) innerIters_ = d;

    this->om_->stream(Debug)
      << " >> stop reason is " << stopReasons_[innerStop_] << endl
      << endl;

  } // end of solveTRSubproblem


#define SIRTR_GET_TEMP_MV(mv,workspace) \
  { \
    TEUCHOS_TEST_FOR_EXCEPTION(workspace.size() == 0,std::logic_error,"SIRTR: Request for workspace could not be honored."); \
    mv = workspace.back(); \
    workspace.pop_back(); \
  }

#define SIRTR_RELEASE_TEMP_MV(mv,workspace) \
  { \
    workspace.push_back(mv); \
    mv = Teuchos::null; \
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Eigensolver iterate() method
  template <class ScalarType, class MV, class OP>
  void SIRTR<ScalarType,MV,OP>::iterate() {

    using Teuchos::RCP;
    using Teuchos::null;
    using Teuchos::tuple;
    using Teuchos::TimeMonitor;
    using std::endl;
    // typedef Teuchos::RCP<const MV> PCMV; // unused
    // typedef Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > PSDM; // unused

    //
    // Allocate/initialize data structures
    //
    if (this->initialized_ == false) {
      this->initialize();
    }

    Teuchos::SerialDenseMatrix<int,ScalarType> AA(this->blockSize_,this->blockSize_),
                                               BB(this->blockSize_,this->blockSize_),
                                               S(this->blockSize_,this->blockSize_);

    // we will often exploit temporarily unused storage for workspace
    // in order to keep it straight and make for clearer code,
    // we will put pointers to available multivectors into the following vector
    // when we need them, we get them out, using a meaningfully-named pointer
    // when we're done, we put them back
    std::vector< RCP<MV> > workspace;
    // we only have 7 multivectors, so that is more than the maximum number that
    // we could use for temp storage
    workspace.reserve(7);

    // set iteration details to invalid, as they don't have any meaning right now
    innerIters_ = -1;
    innerStop_  = UNINITIALIZED;

    // allocate temporary space
    while (this->tester_->checkStatus(this) != Passed) {

      // Print information on current status
      if (this->om_->isVerbosity(Debug)) {
        this->currentStatus( this->om_->stream(Debug) );
      }
      else if (this->om_->isVerbosity(IterationDetails)) {
        this->currentStatus( this->om_->stream(IterationDetails) );
      }

      // increment iteration counter
      this->iter_++;

      // solve the trust-region subproblem
      solveTRSubproblem();
      totalInnerIters_ += innerIters_;

      // perform debugging on eta et al.
      if (this->om_->isVerbosity( Debug ) ) {
        typename RTRBase<ScalarType,MV,OP>::CheckList chk;
        // this is the residual of the model, should still be in the tangent plane
        chk.checkBR  = true;
        chk.checkEta = true;
        this->om_->print( Debug, this->accuracyCheck(chk, "in iterate() after solveTRSubproblem()") );
      }


      //
      // multivectors X, BX (if hasB) and eta contain meaningful information that we need below
      // the others will be sacrificed to temporary storage
      // we are not allowed to reference these anymore, RELEASE_TEMP_MV will clear the pointers
      // the RCP in workspace will keep the MV alive, we will get the MVs back
      // as we need them using GET_TEMP_MV
      //
      // this strategy doesn't cost us much, and it keeps us honest
      //
      TEUCHOS_TEST_FOR_EXCEPTION(workspace.size() != 0,std::logic_error,"SIRTR::iterate(): workspace list should be empty.");
      SIRTR_RELEASE_TEMP_MV(this->delta_ ,workspace);     // workspace size is 1
      SIRTR_RELEASE_TEMP_MV(this->Hdelta_,workspace);     // workspace size is 2
      SIRTR_RELEASE_TEMP_MV(this->R_     ,workspace);     // workspace size is 3
      SIRTR_RELEASE_TEMP_MV(this->Z_     ,workspace);     // workspace size is 4


      // compute the retraction of eta: R_X(eta) = X+eta
      // we must accept, but we will work out of temporary so that we can multiply back into X below
      RCP<MV> XpEta;
      SIRTR_GET_TEMP_MV(XpEta,workspace);                 // workspace size is 3
      {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        TimeMonitor lcltimer( *this->timerLocalUpdate_ );
#endif
        MVT::MvAddMv(ONE,*this->X_,ONE,*this->eta_,*XpEta);
      }

      //
      // perform rayleigh-ritz for XpEta = X+eta
      // save an old copy of f(X) for rho analysis below
      //
      MagnitudeType oldfx = this->fx_;
      int rank, ret;
      rank = this->blockSize_;
      // compute AA = (X+eta)'*A*(X+eta)
      // get temporarily storage for A*(X+eta)
      RCP<MV> AXpEta;
      SIRTR_GET_TEMP_MV(AXpEta,workspace);                // workspace size is 2
      {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        TimeMonitor lcltimer( *this->timerAOp_ );
#endif
        OPT::Apply(*this->AOp_,*XpEta,*AXpEta);
        this->counterAOp_ += this->blockSize_;
      }
      // project A onto X+eta
      {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        TimeMonitor lcltimer( *this->timerLocalProj_ );
#endif
        MVT::MvTransMv(ONE,*XpEta,*AXpEta,AA);
      }
      // compute BB = (X+eta)'*B*(X+eta)
      // get temporary storage for B*(X+eta)
      RCP<MV> BXpEta;
      if (this->hasBOp_) {
        SIRTR_GET_TEMP_MV(BXpEta,workspace);              // workspace size is 1
        {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
          TimeMonitor lcltimer( *this->timerBOp_ );
#endif
          OPT::Apply(*this->BOp_,*XpEta,*BXpEta);
          this->counterBOp_ += this->blockSize_;
        }
        // project B onto X+eta
        {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
          TimeMonitor lcltimer( *this->timerLocalProj_ );
#endif
          MVT::MvTransMv(ONE,*XpEta,*BXpEta,BB);
        }
      }
      else {
        // project I onto X+eta
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        TimeMonitor lcltimer( *this->timerLocalProj_ );
#endif
        MVT::MvTransMv(ONE,*XpEta,*XpEta,BB);
      }
      this->om_->stream(Debug) << "AA: " << std::endl << printMat(AA) << std::endl;;
      this->om_->stream(Debug) << "BB: " << std::endl << printMat(BB) << std::endl;;
      // do the direct solve
      // save old theta first
      std::vector<MagnitudeType> oldtheta(this->theta_);
      {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        TimeMonitor lcltimer( *this->timerDS_ );
#endif
        ret = Utils::directSolver(this->blockSize_,AA,Teuchos::rcpFromRef(BB),S,this->theta_,rank,1);
      }
      this->om_->stream(Debug) << "S: " << std::endl << printMat(S) << std::endl;;
      TEUCHOS_TEST_FOR_EXCEPTION(ret != 0,std::logic_error,"Anasazi::SIRTR::iterate(): failure solving projected eigenproblem after retraction. ret == " << ret << "AA: " << printMat(AA) << std::endl << "BB: " << printMat(BB) << std::endl);
      TEUCHOS_TEST_FOR_EXCEPTION(rank != this->blockSize_,RTRRitzFailure,"Anasazi::SIRTR::iterate(): retracted iterate failed in Ritz analysis. rank == " << rank);

      //
      // order the projected ritz values and vectors
      // this ensures that the ritz vectors produced below are ordered
      {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        TimeMonitor lcltimer( *this->timerSort_ );
#endif
        std::vector<int> order(this->blockSize_);
        // sort the first blockSize_ values in theta_
        this->sm_->sort(this->theta_, Teuchos::rcpFromRef(order), this->blockSize_);   // don't catch exception
        // apply the same ordering to the primitive ritz vectors
        Utils::permuteVectors(order,S);
      }
      //
      // update f(x)
      this->fx_ = std::accumulate(this->theta_.begin(),this->theta_.end(),ZERO);

      //
      // if debugging, do rho analysis before overwriting X,AX,BX
      RCP<MV> AX;
      SIRTR_GET_TEMP_MV(AX,workspace);                   // workspace size is 0
      if (this->om_->isVerbosity( Debug ) ) {
        //
        // compute rho
        //        f(X) - f(X+eta)         f(X) - f(X+eta)
        // rho = ----------------- = -------------------------
        //         m(0) - m(eta)      -<2AX,eta> - .5*<Heta,eta>
        MagnitudeType rhonum, rhoden, mxeta;
        //
        // compute rhonum
        rhonum = oldfx - this->fx_;

        //
        // compute rhoden = -<eta,gradfx> - 0.5 <eta,H*eta>
        //                = -2.0*<eta,AX> - <eta,A*eta> + <eta,B*eta*theta>
        // in three steps        (3)            (1)              (2)
        //
        // first, compute seconder-order decrease in model (steps 1 and 2)
        // get temp storage for second order decrease of model
        //
        // do the first-order decrease last, because we need AX below
        {
          // compute A*eta and then <eta,A*eta>
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
          TimeMonitor lcltimer( *this->timerAOp_ );
#endif
          OPT::Apply(*this->AOp_,*this->eta_,*AX);
          this->counterAOp_ += this->blockSize_;
        }
        // compute A part of second order decrease into rhoden
        rhoden = -RTRBase<ScalarType,MV,OP>::ginner(*this->eta_,*AX);
        if (this->hasBOp_) {
          // compute B*eta into AX
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
          TimeMonitor lcltimer( *this->timerBOp_ );
#endif
          OPT::Apply(*this->BOp_,*this->eta_,*AX);
          this->counterBOp_ += this->blockSize_;
        }
        else {
          // put B*eta==eta into AX
          MVT::MvAddMv(ONE,*this->eta_,ZERO,*this->eta_,*AX);
        }
        // we need this below for computing individual rho, get it now
        std::vector<MagnitudeType> eBe(this->blockSize_);
        RTRBase<ScalarType,MV,OP>::ginnersep(*this->eta_,*AX,eBe);
        // scale B*eta by theta
        {
          std::vector<ScalarType> oldtheta_complex(oldtheta.begin(),oldtheta.end());
          MVT::MvScale( *AX, oldtheta_complex);
        }
        // accumulate B part of second order decrease into rhoden
        rhoden += RTRBase<ScalarType,MV,OP>::ginner(*this->eta_,*AX);
        {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
          TimeMonitor lcltimer( *this->timerAOp_ );
#endif
          OPT::Apply(*this->AOp_,*this->X_,*AX);
          this->counterAOp_ += this->blockSize_;
        }
        // accumulate first-order decrease of model into rhoden
        rhoden += -2.0*RTRBase<ScalarType,MV,OP>::ginner(*AX,*this->eta_);

        mxeta = oldfx - rhoden;
        this->rho_ = rhonum / rhoden;
        this->om_->stream(Debug)
          << " >> old f(x) is : " << oldfx << endl
          << " >> new f(x) is : " << this->fx_ << endl
          << " >> m_x(eta) is : " << mxeta << endl
          << " >> rhonum is   : " << rhonum << endl
          << " >> rhoden is   : " << rhoden << endl
          << " >> rho is      : " << this->rho_ << endl;
        // compute individual rho
        for (int j=0; j<this->blockSize_; ++j) {
          this->om_->stream(Debug)
            << " >> rho[" << j << "]     : " << 1.0/(1.0+eBe[j]) << endl;
        }
      }

      // compute Ritz vectors back into X,BX,AX
      {
        // release const views to X, BX
        this->X_  = Teuchos::null;
        this->BX_ = Teuchos::null;
        // get non-const views
        std::vector<int> ind(this->blockSize_);
        for (int i=0; i<this->blockSize_; ++i) ind[i] = this->numAuxVecs_+i;
        Teuchos::RCP<MV> X, BX;
        X = MVT::CloneViewNonConst(*this->V_,ind);
        if (this->hasBOp_) {
          BX = MVT::CloneViewNonConst(*this->BV_,ind);
        }
        // compute ritz vectors, A,B products into X,AX,BX
        {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
          TimeMonitor lcltimer( *this->timerLocalUpdate_ );
#endif
          MVT::MvTimesMatAddMv(ONE,* XpEta,S,ZERO,*X);
          MVT::MvTimesMatAddMv(ONE,*AXpEta,S,ZERO,*AX);
          if (this->hasBOp_) {
            MVT::MvTimesMatAddMv(ONE,*BXpEta,S,ZERO,*BX);
          }
        }
        // clear non-const views, restore const views
        X  = Teuchos::null;
        BX = Teuchos::null;
        this->X_  = MVT::CloneView(static_cast<const MV&>(*this->V_ ),ind);
        this->BX_ = MVT::CloneView(static_cast<const MV&>(*this->BV_),ind);
      }
      //
      // return XpEta and BXpEta to temp storage
      SIRTR_RELEASE_TEMP_MV(XpEta,workspace);             // workspace size is 1
      SIRTR_RELEASE_TEMP_MV(AXpEta,workspace);            // workspace size is 2
      if (this->hasBOp_) {
        SIRTR_RELEASE_TEMP_MV(BXpEta,workspace);          // workspace size is 3
      }

      //
      // solveTRSubproblem destroyed R, we must recompute it
      // compute R = AX - BX*theta
      //
      // get R back from temp storage
      SIRTR_GET_TEMP_MV(this->R_,workspace);              // workspace size is 2
      {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        TimeMonitor lcltimer( *this->timerCompRes_ );
#endif
        MVT::MvAddMv( ONE, *this->BX_, ZERO, *this->BX_, *this->R_ );
        {
          std::vector<ScalarType> theta_comp(this->theta_.begin(),this->theta_.end());
          MVT::MvScale( *this->R_, theta_comp );
        }
        MVT::MvAddMv( ONE, *AX, -ONE, *this->R_, *this->R_ );
      }
      //
      // R has been updated; mark the norms as out-of-date
      this->Rnorms_current_ = false;
      this->R2norms_current_ = false;

      //
      // we are done with AX, release it
      SIRTR_RELEASE_TEMP_MV(AX,workspace);                // workspace size is 3
      //
      // get data back for delta, Hdelta and Z
      SIRTR_GET_TEMP_MV(this->delta_,workspace);          // workspace size is 2
      SIRTR_GET_TEMP_MV(this->Hdelta_,workspace);         // workspace size is 1
      SIRTR_GET_TEMP_MV(this->Z_,workspace);              // workspace size is 0

      //
      // When required, monitor some orthogonalities
      if (this->om_->isVerbosity( Debug ) ) {
        // Check almost everything here
        typename RTRBase<ScalarType,MV,OP>::CheckList chk;
        chk.checkX = true;
        chk.checkBX = true;
        chk.checkR = true;
        this->om_->print( Debug, this->accuracyCheck(chk, "after local update") );
      }
      else if (this->om_->isVerbosity( OrthoDetails )) {
        typename RTRBase<ScalarType,MV,OP>::CheckList chk;
        chk.checkX = true;
        chk.checkR = true;
        this->om_->print( OrthoDetails, this->accuracyCheck(chk, "after local update") );
      }

    } // end while (statusTest == false)

  } // end of iterate()


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Print the current status of the solver
  template <class ScalarType, class MV, class OP>
  void
  SIRTR<ScalarType,MV,OP>::currentStatus(std::ostream &os)
  {
    using std::endl;
    os.setf(std::ios::scientific, std::ios::floatfield);
    os.precision(6);
    os <<endl;
    os <<"================================================================================" << endl;
    os << endl;
    os <<"                         SIRTR Solver Status" << endl;
    os << endl;
    os <<"The solver is "<<(this->initialized_ ? "initialized." : "not initialized.") << endl;
    os <<"The number of iterations performed is " << this->iter_       << endl;
    os <<"The current block size is             " << this->blockSize_  << endl;
    os <<"The number of auxiliary vectors is    " << this->numAuxVecs_ << endl;
    os <<"The number of operations A*x    is " << this->counterAOp_   << endl;
    os <<"The number of operations B*x    is " << this->counterBOp_    << endl;
    os <<"The number of operations B*x by the orthomanager is " << this->orthman_->getOpCounter() << endl;
    os <<"The number of operations Prec*x is " << this->counterPrec_ << endl;
    os <<"Parameter rho_prime is  " << rho_prime_ << endl;
    os <<"Inner stopping condition was " << stopReasons_[innerStop_] << endl;
    os <<"Number of inner iterations was " << innerIters_ << endl;
    os <<"Total number of inner iterations is " << totalInnerIters_ << endl;
    os <<"f(x) is " << this->fx_ << endl;

    os.setf(std::ios_base::right, std::ios_base::adjustfield);

    if (this->initialized_) {
      os << endl;
      os <<"CURRENT EIGENVALUE ESTIMATES             "<<endl;
      os << std::setw(20) << "Eigenvalue"
         << std::setw(20) << "Residual(B)"
         << std::setw(20) << "Residual(2)"
         << endl;
      os <<"--------------------------------------------------------------------------------"<<endl;
      for (int i=0; i<this->blockSize_; i++) {
        os << std::setw(20) << this->theta_[i];
        if (this->Rnorms_current_) os << std::setw(20) << this->Rnorms_[i];
        else os << std::setw(20) << "not current";
        if (this->R2norms_current_) os << std::setw(20) << this->R2norms_[i];
        else os << std::setw(20) << "not current";
        os << endl;
      }
    }
    os <<"================================================================================" << endl;
    os << endl;
  }


} // end Anasazi namespace

#endif // ANASAZI_SIRTR_HPP
