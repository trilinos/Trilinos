// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//! \file AnasaziIRTR.hpp


#ifndef ANASAZI_IRTR_HPP
#define ANASAZI_IRTR_HPP

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

/*!     \class Anasazi::IRTR

        IRTR is a caching implementation of the Implicit Riemannian Trust-Region
        (%IRTR) eigensolver.

        The solver uses between 10 and 13 blocks of vectors, compared to the
        requirements by SIRTR of 6 to 8 blocks of vectors. The base requirement
        is 10 blocks of vectors, where a block of vectors contains a number of vectors equal to the
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
  class IRTR : public RTRBase<ScalarType,MV,OP> {
  public:

    //! @name Constructor/Destructor
    //@{

    /*! \brief %IRTR constructor with eigenproblem, solver utilities, and parameter list of solver options.
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
    IRTR( const Teuchos::RCP<Eigenproblem<ScalarType,MV,OP> >    &problem,
          const Teuchos::RCP<SortManager<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> >           &sorter,
          const Teuchos::RCP<OutputManager<ScalarType> >         &printer,
          const Teuchos::RCP<StatusTest<ScalarType,MV,OP> >      &tester,
          const Teuchos::RCP<GenOrthoManager<ScalarType,MV,OP> > &ortho,
          Teuchos::ParameterList                                 &params
        );

    //! %IRTR destructor
    virtual ~IRTR() {};

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
    //
    // use 2D subspace acceleration of X+Eta to generate new iterate?
    bool useSA_;
  };




  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor
  template <class ScalarType, class MV, class OP>
  IRTR<ScalarType,MV,OP>::IRTR(
        const Teuchos::RCP<Eigenproblem<ScalarType,MV,OP> >    &problem,
        const Teuchos::RCP<SortManager<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> >           &sorter,
        const Teuchos::RCP<OutputManager<ScalarType> >         &printer,
        const Teuchos::RCP<StatusTest<ScalarType,MV,OP> >      &tester,
        const Teuchos::RCP<GenOrthoManager<ScalarType,MV,OP> > &ortho,
        Teuchos::ParameterList                                 &params
        ) :
    RTRBase<ScalarType,MV,OP>(problem,sorter,printer,tester,ortho,params,"IRTR",false),
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
                       "Anasazi::IRTR::constructor: rho_prime must be in (0,1).");

    useSA_ = params.get<bool>("Use SA",false);
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // TR subproblem solver
  //
  // FINISH:
  //   define pre- and post-conditions
  //
  // POST:
  //   delta_,Adelta_,Bdelta_,Hdelta_ undefined
  //
  template <class ScalarType, class MV, class OP>
  void IRTR<ScalarType,MV,OP>::solveTRSubproblem() {

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
    MVT::MvInit(*this->Aeta_,0.0);
    if (this->hasBOp_) {
      MVT::MvInit(*this->Beta_,0.0);
    }

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
      RCP<MV> Heta = MVT::Clone(*this->eta_,this->blockSize_);
      // put 2*A*e - 2*B*e*theta --> He
      MVT::MvAddMv(ONE,*this->Beta_,ZERO,*this->Beta_,*Heta);
      std::vector<ScalarType> theta_comp(this->theta_.begin(),this->theta_.end());
      MVT::MvScale(*Heta,theta_comp);
      if (this->leftMost_) {
        MVT::MvAddMv( 2.0,*this->Aeta_,-2.0,*Heta,*Heta); // Heta = 2*Aeta + (-2)*Beta*theta
      }
      else {
        MVT::MvAddMv(-2.0,*this->Aeta_, 2.0,*Heta,*Heta); // Heta = (-2)*Aeta + 2*Beta*theta
      }
      // apply projector
      {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        TimeMonitor lcltimer( *this->timerOrtho_ );
#endif
        this->orthman_->projectGen(
            *Heta,                                                 // operating on Heta
            tuple<PCMV>(this->BV_),tuple<PCMV>(this->V_),false,    // P_{BV,V}, and <BV,V>_B != I
            tuple<PSDM>(null),                                     // don't care about coeffs
            null,tuple<PCMV>(null), tuple<PCMV>(this->BV_));       // don't have B*BV, but do have B*V
      }

      std::vector<MagnitudeType> eg(this->blockSize_),
                                eHe(this->blockSize_);
      RTRBase<ScalarType,MV,OP>::ginnersep(*this->eta_,*this->AX_,eg);
      RTRBase<ScalarType,MV,OP>::ginnersep(*this->eta_,*Heta,eHe);
      if (this->leftMost_) {
        for (int j=0; j<this->blockSize_; ++j) {
          eg[j] = this->theta_[j] + 2*eg[j] + .5*eHe[j];
        }
      }
      else {
        for (int j=0; j<this->blockSize_; ++j) {
          eg[j] = -this->theta_[j] - 2*eg[j] + .5*eHe[j];
        }
      }
      this->om_->stream(Debug)
        << " Debugging checks: IRTR inner iteration " << innerIters_ << endl
        << " >> m_X(eta) : " << std::accumulate(eg.begin(),eg.end(),0.0) << endl;
      for (int j=0; j<this->blockSize_; ++j) {
        this->om_->stream(Debug)
          << " >> m_X(eta_" << j << ") : " << eg[j] << endl;
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
        OPT::Apply(*this->AOp_,*this->delta_,*this->Adelta_);
        this->counterAOp_ += this->blockSize_;
      }
      if (this->hasBOp_) {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        TimeMonitor lcltimer( *this->timerBOp_ );
#endif
        OPT::Apply(*this->BOp_,*this->delta_,*this->Bdelta_);
        this->counterBOp_ += this->blockSize_;
      }
      else {
        MVT::MvAddMv(ONE,*this->delta_,ZERO,*this->delta_,*this->Bdelta_);
      }
      // compute <eta,B*delta> and <delta,B*delta>
      // these will be needed below
      RTRBase<ScalarType,MV,OP>::ginnersep(*this->eta_  ,*this->Bdelta_,eBd);
      RTRBase<ScalarType,MV,OP>::ginnersep(*this->delta_,*this->Bdelta_,dBd);
      // put 2*A*d - 2*B*d*theta --> Hd
      MVT::MvAddMv(ONE,*this->Bdelta_,ZERO,*this->Bdelta_,*this->Hdelta_);
      {
        std::vector<ScalarType> theta_comp(this->theta_.begin(),this->theta_.end());
        MVT::MvScale(*this->Hdelta_,theta_comp);
      }
      if (this->leftMost_) {
        MVT::MvAddMv( 2.0,*this->Adelta_,-2.0,*this->Hdelta_,*this->Hdelta_);
      }
      else {
        MVT::MvAddMv(-2.0,*this->Adelta_, 2.0,*this->Hdelta_,*this->Hdelta_);
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
        MVT::MvScale(*this->Adelta_,alpha_comp);
        if (this->hasBOp_) {
          MVT::MvScale(*this->Bdelta_,alpha_comp);
        }
        MVT::MvScale(*this->Hdelta_,alpha_comp);
      }
      MVT::MvAddMv(ONE,*this->delta_ ,ONE,*this->eta_ ,*this->eta_);
      MVT::MvAddMv(ONE,*this->Adelta_,ONE,*this->Aeta_,*this->Aeta_);
      if (this->hasBOp_) {
        MVT::MvAddMv(ONE,*this->Bdelta_,ONE,*this->Beta_,*this->Beta_);
      }

      // store new eBe
      eBe = new_eBe;

      //
      // print some debugging info
      if (this->om_->isVerbosity(Debug)) {
        RCP<MV> Heta = MVT::Clone(*this->eta_,this->blockSize_);
        // put 2*A*e - 2*B*e*theta --> He
        MVT::MvAddMv(ONE,*this->Beta_,ZERO,*this->Beta_,*Heta);
        {
          std::vector<ScalarType> theta_comp(this->theta_.begin(),this->theta_.end());
          MVT::MvScale(*Heta,theta_comp);
        }
        if (this->leftMost_) {
          MVT::MvAddMv( 2.0,*this->Aeta_,-2.0,*Heta,*Heta);
        }
        else {
          MVT::MvAddMv(-2.0,*this->Aeta_, 2.0,*Heta,*Heta);
        }
        // apply projector
        {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
          TimeMonitor lcltimer( *this->timerOrtho_ );
#endif
          this->orthman_->projectGen(
              *Heta,                                                // operating on Heta
              tuple<PCMV>(this->BV_),tuple<PCMV>(this->V_),false,   // P_{BV,V}, and <BV,V>_B != I
              tuple<PSDM>(null),                                    // don't care about coeffs
              null,tuple<PCMV>(null), tuple<PCMV>(this->BV_));      // don't have B*BV, but do have B*V
        }

        std::vector<MagnitudeType> eg(this->blockSize_),
                                   eHe(this->blockSize_);
        RTRBase<ScalarType,MV,OP>::ginnersep(*this->eta_,*this->AX_,eg);
        RTRBase<ScalarType,MV,OP>::ginnersep(*this->eta_,*Heta,eHe);
        if (this->leftMost_) {
          for (int j=0; j<this->blockSize_; ++j) {
            eg[j] = this->theta_[j] + 2*eg[j] + .5*eHe[j];
          }
        }
        else {
          for (int j=0; j<this->blockSize_; ++j) {
            eg[j] = -this->theta_[j] - 2*eg[j] + .5*eHe[j];
          }
        }
        this->om_->stream(Debug)
          << " Debugging checks: IRTR inner iteration " << innerIters_ << endl
          << " >> m_X(eta) : " << std::accumulate(eg.begin(),eg.end(),0.0) << endl;
        for (int j=0; j<this->blockSize_; ++j) {
          this->om_->stream(Debug)
            << " >> m_X(eta_" << j << ") : " << eg[j] << endl;
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


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Eigensolver iterate() method
  template <class ScalarType, class MV, class OP>
  void IRTR<ScalarType,MV,OP>::iterate() {

    using Teuchos::RCP;
    using Teuchos::null;
    using Teuchos::tuple;
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
    using Teuchos::TimeMonitor;
#endif
    using std::endl;
    //typedef Teuchos::RCP<const MV> PCMV; // unused
    //typedef Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > PSDM; // unused

    //
    // Allocate/initialize data structures
    //
    if (this->initialized_ == false) {
      this->initialize();
    }

    Teuchos::SerialDenseMatrix<int,ScalarType> AA, BB, S;
    if (useSA_ == true) {
      AA.reshape(2*this->blockSize_,2*this->blockSize_);
      BB.reshape(2*this->blockSize_,2*this->blockSize_);
      S.reshape(2*this->blockSize_,2*this->blockSize_);
    }
    else {
      AA.reshape(this->blockSize_,this->blockSize_);
      BB.reshape(this->blockSize_,this->blockSize_);
      S.reshape(this->blockSize_,this->blockSize_);
    }

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
        chk.checkAEta = true;
        chk.checkBEta = true;
        this->om_->print( Debug, this->accuracyCheck(chk, "in iterate() after solveTRSubproblem()") );
        this->om_->stream(Debug)
          << " >> norm(Eta) : " << MAT::squareroot(RTRBase<ScalarType,MV,OP>::ginner(*this->eta_)) << endl
          << endl;
      }

      if (useSA_ == false) {
        // compute the retraction of eta: R_X(eta) = X+eta
        // we must accept, but we will work out of delta so that we can multiply back into X below
        {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
          TimeMonitor lcltimer( *this->timerLocalUpdate_ );
#endif
          MVT::MvAddMv(ONE,*this->X_,ONE,*this->eta_,*this->delta_);
          MVT::MvAddMv(ONE,*this->AX_,ONE,*this->Aeta_,*this->Adelta_);
          if (this->hasBOp_) {
            MVT::MvAddMv(ONE,*this->BX_,ONE,*this->Beta_,*this->Bdelta_);
          }
        }
        // perform some debugging on X+eta
        if (this->om_->isVerbosity( Debug ) ) {
          // X^T B X = I
          // X^T B Eta = 0
          // (X+Eta)^T B (X+Eta) = I + Eta^T B Eta
          Teuchos::SerialDenseMatrix<int,ScalarType> XE(this->blockSize_,this->blockSize_),
            E(this->blockSize_,this->blockSize_);
          MVT::MvTransMv(ONE,*this->delta_,*this->Bdelta_,XE);
          MVT::MvTransMv(ONE,*this->eta_,*this->Beta_,E);
          this->om_->stream(Debug)
            << " >> Error in AX+AEta == A(X+Eta) : " << Utils::errorEquality(*this->delta_,*this->Adelta_,this->AOp_) << endl
            << " >> Error in BX+BEta == B(X+Eta) : " << Utils::errorEquality(*this->delta_,*this->Bdelta_,this->BOp_) << endl
            << " >> norm( (X+Eta)^T B (X+Eta) )  : " << XE.normFrobenius() << endl
            << " >> norm( Eta^T B Eta )          : " << E.normFrobenius() << endl
            << endl;
        }
      }

      //
      // perform rayleigh-ritz for X+eta or [X,eta] according to useSA_
      // save an old copy of f(X) for rho analysis below
      //
      MagnitudeType oldfx = this->fx_;
      std::vector<MagnitudeType> oldtheta(this->theta_), newtheta(AA.numRows());
      int rank, ret;
      rank = AA.numRows();
      if (useSA_ == true) {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        TimeMonitor lcltimer( *this->timerLocalProj_ );
#endif
        Teuchos::SerialDenseMatrix<int,ScalarType> AA11(Teuchos::View,AA,this->blockSize_,this->blockSize_,0,0),
                                                   AA12(Teuchos::View,AA,this->blockSize_,this->blockSize_,0,this->blockSize_),
                                                   AA22(Teuchos::View,AA,this->blockSize_,this->blockSize_,this->blockSize_,this->blockSize_);
        Teuchos::SerialDenseMatrix<int,ScalarType> BB11(Teuchos::View,BB,this->blockSize_,this->blockSize_,0,0),
                                                   BB12(Teuchos::View,BB,this->blockSize_,this->blockSize_,0,this->blockSize_),
                                                   BB22(Teuchos::View,BB,this->blockSize_,this->blockSize_,this->blockSize_,this->blockSize_);
        MVT::MvTransMv(ONE,*this->X_  ,*this->AX_  ,AA11);
        MVT::MvTransMv(ONE,*this->X_  ,*this->Aeta_,AA12);
        MVT::MvTransMv(ONE,*this->eta_,*this->Aeta_,AA22);
        MVT::MvTransMv(ONE,*this->X_  ,*this->BX_  ,BB11);
        MVT::MvTransMv(ONE,*this->X_  ,*this->Beta_,BB12);
        MVT::MvTransMv(ONE,*this->eta_,*this->Beta_,BB22);
      }
      else {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        TimeMonitor lcltimer( *this->timerLocalProj_ );
#endif
        MVT::MvTransMv(ONE,*this->delta_,*this->Adelta_,AA);
        MVT::MvTransMv(ONE,*this->delta_,*this->Bdelta_,BB);
      }
      this->om_->stream(Debug) << "AA: " << std::endl << printMat(AA) << std::endl;;
      this->om_->stream(Debug) << "BB: " << std::endl << printMat(BB) << std::endl;;
      {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        TimeMonitor lcltimer( *this->timerDS_ );
#endif
        ret = Utils::directSolver(AA.numRows(),AA,Teuchos::rcpFromRef(BB),S,newtheta,rank,1);
      }
      this->om_->stream(Debug) << "S: " << std::endl << printMat(S) << std::endl;;
      TEUCHOS_TEST_FOR_EXCEPTION(ret != 0,std::logic_error,"Anasazi::IRTR::iterate(): failure solving projected eigenproblem after retraction. ret == " << ret);
      TEUCHOS_TEST_FOR_EXCEPTION(rank != AA.numRows(),RTRRitzFailure,"Anasazi::IRTR::iterate(): retracted iterate failed in Ritz analysis. rank == " << rank);

      //
      // order the projected ritz values and vectors
      // this ensures that the ritz vectors produced below are ordered
      {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        TimeMonitor lcltimer( *this->timerSort_ );
#endif
        std::vector<int> order(newtheta.size());
        // sort the values in newtheta
        this->sm_->sort(newtheta, Teuchos::rcpFromRef(order), -1);   // don't catch exception
        // apply the same ordering to the primitive ritz vectors
        Utils::permuteVectors(order,S);
      }
      //
      // save the first blockSize values into this->theta_
      std::copy(newtheta.begin(), newtheta.begin()+this->blockSize_, this->theta_.begin());
      //
      // update f(x)
      this->fx_ = std::accumulate(this->theta_.begin(),this->theta_.end(),ZERO);

      //
      // if debugging, do rho analysis before overwriting X,AX,BX
      if (this->om_->isVerbosity( Debug ) ) {
        //
        // compute rho
        //        f(X) - f(newX)           f(X) - f(newX)
        // rho = ---------------- = ---------------------------
        //         m(0) - m(eta)    -<2AX,eta> - .5*<Heta,eta>
        //
        //            f(X) - f(newX)
        //     = ---------------------------------------
        //        -<2AX,eta> - <eta,Aeta> + <eta,Beta XAX>
        //
        MagnitudeType rhonum, rhoden, mxeta;
        std::vector<MagnitudeType> eBe(this->blockSize_);
        RTRBase<ScalarType,MV,OP>::ginnersep(*this->eta_,*this->Beta_,eBe);
        //
        // compute rhonum
        rhonum = oldfx - this->fx_;
        //
        // compute rhoden
        rhoden = -2.0*RTRBase<ScalarType,MV,OP>::ginner(*this->AX_  ,*this->eta_)
                 -RTRBase<ScalarType,MV,OP>::ginner(*this->Aeta_,*this->eta_);
        for (int i=0; i<this->blockSize_; ++i) {
          rhoden += eBe[i]*oldtheta[i];
        }
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

      // form new X as Ritz vectors, using the primitive Ritz vectors in S and
      // either [X eta] or X+eta
      // we will clear the const views of X,BX into V,BV and
      // work from non-const temporary views
      {
        // release const views to X, BX
        // get non-const views
        this->X_  = Teuchos::null;
        this->BX_ = Teuchos::null;
        std::vector<int> ind(this->blockSize_);
        for (int i=0; i<this->blockSize_; ++i) ind[i] = this->numAuxVecs_+i;
        Teuchos::RCP<MV> X, BX;
        X = MVT::CloneViewNonConst(*this->V_,ind);
        if (this->hasBOp_) {
          BX = MVT::CloneViewNonConst(*this->BV_,ind);
        }
        if (useSA_ == false) {
          // multiply delta=(X+eta),Adelta=...,Bdelta=...
          // by primitive Ritz vectors back into X,AX,BX
          // compute ritz vectors, A,B products into X,AX,BX
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
          TimeMonitor lcltimer( *this->timerLocalUpdate_ );
#endif
          MVT::MvTimesMatAddMv(ONE,*this->delta_,S,ZERO,*X);
          MVT::MvTimesMatAddMv(ONE,*this->Adelta_,S,ZERO,*this->AX_);
          if (this->hasBOp_) {
            MVT::MvTimesMatAddMv(ONE,*this->Bdelta_,S,ZERO,*BX);
          }
        }
        else {
          // compute ritz vectors, A,B products into X,AX,BX
          // currently, X in X and eta in eta
          // compute each result into delta, then copy to appropriate place
          // decompose S into [Sx;Se]
          Teuchos::SerialDenseMatrix<int,ScalarType> Sx(Teuchos::View,S,this->blockSize_,this->blockSize_,0,0),
                                                     Se(Teuchos::View,S,this->blockSize_,this->blockSize_,this->blockSize_,0);
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
          TimeMonitor lcltimer( *this->timerLocalUpdate_ );
#endif
          // X = [X eta] S = X*Sx + eta*Se
          MVT::MvTimesMatAddMv(ONE,*      X   ,Sx,ZERO,*this->delta_);
          MVT::MvTimesMatAddMv(ONE,*this->eta_,Se,ONE ,*this->delta_);
          MVT::MvAddMv(ONE,*this->delta_,ZERO,*this->delta_,*X);
          // AX = [AX Aeta] S = AX*Sx + Aeta*Se
          MVT::MvTimesMatAddMv(ONE,*this->AX_  ,Sx,ZERO,*this->delta_);
          MVT::MvTimesMatAddMv(ONE,*this->Aeta_,Se,ONE ,*this->delta_);
          MVT::MvAddMv(ONE,*this->delta_,ZERO,*this->delta_,*this->AX_);
          if (this->hasBOp_) {
            // BX = [BX Beta] S = BX*Sx + Beta*Se
            MVT::MvTimesMatAddMv(ONE,*      BX   ,Sx,ZERO,*this->delta_);
            MVT::MvTimesMatAddMv(ONE,*this->Beta_,Se,ONE ,*this->delta_);
            MVT::MvAddMv(ONE,*this->delta_,ZERO,*this->delta_,*BX);
          }
        }
        // clear non-const views, restore const views
        X  = Teuchos::null;
        BX = Teuchos::null;
        this->X_  = MVT::CloneView(static_cast<const MV&>(*this->V_ ),ind);
        this->BX_ = MVT::CloneView(static_cast<const MV&>(*this->BV_),ind);
      }

      //
      // update residual, R = AX - BX*theta
      {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        TimeMonitor lcltimer( *this->timerCompRes_ );
#endif
        MVT::MvAddMv( ONE, *this->BX_, ZERO, *this->BX_, *this->R_ );
        std::vector<ScalarType> theta_comp(this->theta_.begin(),this->theta_.end());
        MVT::MvScale( *this->R_, theta_comp );
        MVT::MvAddMv( ONE, *this->AX_, -ONE, *this->R_, *this->R_ );
      }
      //
      // R has been updated; mark the norms as out-of-date
      this->Rnorms_current_ = false;
      this->R2norms_current_ = false;


      //
      // When required, monitor some orthogonalities
      if (this->om_->isVerbosity( Debug ) ) {
        // Check almost everything here
        typename RTRBase<ScalarType,MV,OP>::CheckList chk;
        chk.checkX = true;
        chk.checkAX = true;
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
  IRTR<ScalarType,MV,OP>::currentStatus(std::ostream &os)
  {
    using std::endl;
    os.setf(std::ios::scientific, std::ios::floatfield);
    os.precision(6);
    os <<endl;
    os <<"================================================================================" << endl;
    os << endl;
    os <<"                             IRTR Solver Status" << endl;
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

#endif // ANASAZI_IRTR_HPP
