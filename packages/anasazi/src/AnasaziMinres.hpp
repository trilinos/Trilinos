// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file AnasaziMinres.hpp
 *  An implementation of pseudo-block MINRES capable of working with different operators for each RHS.
*/

#ifndef ANASAZI_MINRES_HPP
#define ANASAZI_MINRES_HPP

#include "AnasaziConfigDefs.hpp"
#include "Teuchos_TimeMonitor.hpp"

namespace Anasazi {
namespace Experimental {

template <class Scalar, class MV, class OP>
class PseudoBlockMinres
{
  typedef Anasazi::MultiVecTraits<Scalar,MV>         MVT;
  typedef Anasazi::OperatorTraits<Scalar,MV,OP>      OPT;
  const Scalar ONE; 
  const Scalar ZERO;

public:
  // Constructor
  PseudoBlockMinres(Teuchos::RCP<OP> A, Teuchos::RCP<OP> Prec = Teuchos::null);

  // Set tolerance and maximum iterations
  void setTol(const std::vector<Scalar>& tolerances) { tolerances_ = tolerances; };
  void setMaxIter(const int maxIt) { maxIt_ = maxIt; };

  // Set solution and RHS
  void setSol(Teuchos::RCP<MV> X) { X_ = X; };
  void setRHS(Teuchos::RCP<const MV> B) { B_ = B; };

  // Solve the linear system
  void solve();

private:
  Teuchos::RCP<OP> A_;
  Teuchos::RCP<OP> Prec_;
  Teuchos::RCP<MV> X_;
  Teuchos::RCP<const MV> B_;
  std::vector<Scalar> tolerances_;
  int maxIt_;

  void symOrtho(Scalar a, Scalar b, Scalar *c, Scalar *s, Scalar *r);

#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
  Teuchos::RCP<Teuchos::Time> AddTime_, ApplyOpTime_, ApplyPrecTime_, AssignTime_, DotTime_, LockTime_, NormTime_, ScaleTime_, TotalTime_;
#endif
};



template <class Scalar, class MV, class OP>
PseudoBlockMinres<Scalar,MV,OP>::PseudoBlockMinres(Teuchos::RCP<OP> A, Teuchos::RCP<OP> Prec) : 
  ONE(Teuchos::ScalarTraits<Scalar>::one()), 
  ZERO(Teuchos::ScalarTraits<Scalar>::zero()) 
{
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
  AddTime_ = Teuchos::TimeMonitor::getNewTimer("Anasazi: PseudoBlockMinres::Add Multivectors"); 
  ApplyOpTime_ = Teuchos::TimeMonitor::getNewTimer("Anasazi: PseudoBlockMinres::Apply Operator");
  ApplyPrecTime_ = Teuchos::TimeMonitor::getNewTimer("Anasazi: PseudoBlockMinres::Apply Preconditioner");
  AssignTime_ = Teuchos::TimeMonitor::getNewTimer("Anasazi: PseudoBlockMinres::Assignment (no locking)");
  DotTime_ = Teuchos::TimeMonitor::getNewTimer("Anasazi: PseudoBlockMinres::Compute Dot Product");
  LockTime_ = Teuchos::TimeMonitor::getNewTimer("Anasazi: PseudoBlockMinres::Lock Converged Vectors");
  NormTime_ = Teuchos::TimeMonitor::getNewTimer("Anasazi: PseudoBlockMinres::Compute Norm");
  ScaleTime_ = Teuchos::TimeMonitor::getNewTimer("Anasazi: PseudoBlockMinres::Scale Multivector");
  TotalTime_ = Teuchos::TimeMonitor::getNewTimer("Anasazi: PseudoBlockMinres::Total Time");
#endif

  A_ = A;
  Prec_ = Prec;
  maxIt_ = 0;
}



template <class Scalar, class MV, class OP>
void PseudoBlockMinres<Scalar,MV,OP>::solve()
{
  #ifdef ANASAZI_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor outertimer( *TotalTime_ );
  #endif

  // Get number of vectors
  int ncols = MVT::GetNumberVecs(*B_);
  int newNumConverged;

  // Declare some variables
  std::vector<Scalar> alpha(ncols), beta, beta1(ncols), phibar, oldBeta(ncols,ZERO), epsln(ncols,ZERO), cs(ncols,-ONE), sn(ncols,ZERO), dbar(ncols,ZERO), oldeps(ncols), delta(ncols), gbar(ncols), phi(ncols), gamma(ncols), tmpvec(ncols);
  std::vector<int> indicesToRemove, newlyConverged, unconvergedIndices(ncols);
  Teuchos::RCP<MV> V, Y, R1, R2, W, W1, W2, locX, scaleHelper, swapHelper;

  // Allocate space for multivectors
  V = MVT::Clone(*B_, ncols);
  Y = MVT::Clone(*B_, ncols);
  R1 = MVT::Clone(*B_, ncols);
  R2 = MVT::Clone(*B_, ncols);
  W = MVT::Clone(*B_, ncols);
  W1 = MVT::Clone(*B_, ncols);
  W2 = MVT::Clone(*B_, ncols);
  scaleHelper = MVT::Clone(*B_, ncols);

  // Reserve space for arrays
  indicesToRemove.reserve(ncols);
  newlyConverged.reserve(ncols);

  // Initialize array of unconverged indices
  for(int i=0; i<ncols; i++)
    unconvergedIndices[i] = i;

  // Get a local copy of X
  // We want the vectors to remain contiguous even as things converge
  locX = MVT::CloneCopy(*X_);

  // Initialize residuals
  // R1 = B - AX
  {
    #ifdef ANASAZI_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor lcltimer( *ApplyOpTime_ );
    #endif
    OPT::Apply(*A_,*locX,*R1);
  }
  {
    #ifdef ANASAZI_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor lcltimer( *AddTime_ );
    #endif
    MVT::MvAddMv(ONE, *B_, -ONE, *R1, *R1);
  }

  // R2 = R1
  {
    #ifdef ANASAZI_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor lcltimer( *AssignTime_ );
    #endif
    MVT::Assign(*R1,*R2);
  }

  // Initialize the W's to 0.
  MVT::MvInit (*W);
  MVT::MvInit (*W2);

  // Y = M\R1 (preconditioned residual)
  if(Prec_ != Teuchos::null)
  {
    #ifdef ANASAZI_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor lcltimer( *ApplyPrecTime_ );
    #endif

    OPT::Apply(*Prec_,*R1,*Y);
  }
  else
  {
    #ifdef ANASAZI_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor lcltimer( *AssignTime_ );
    #endif
    MVT::Assign(*R1,*Y);
  }

  // beta1 = sqrt(Y'*R1)
  {
    #ifdef ANASAZI_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor lcltimer( *DotTime_ );
    #endif
    MVT::MvDot(*Y,*R1, beta1);

    for(size_t i=0; i<beta1.size(); i++)
      beta1[i] = sqrt(beta1[i]);
  }

  // beta = beta1
  beta = beta1;

  // phibar = beta1
  phibar = beta1;

  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Begin Lanczos iterations
  for(int iter=1; iter <= maxIt_; iter++)
  {
    // Test convergence
    indicesToRemove.clear();
    for(int i=0; i<ncols; i++)
    {
      Scalar relres = phibar[i]/beta1[i];
//      std::cout << "relative residual[" << unconvergedIndices[i] << "] at iteration " << iter << " = " << relres << std::endl;

      // If the vector converged, mark it for termination
      // Make sure to add it to 
      if(relres < tolerances_[i])
      {
        indicesToRemove.push_back(i);
      }
    }

    // Check whether anything converged
    newNumConverged = indicesToRemove.size();
    if(newNumConverged > 0)
    {
      #ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor lcltimer( *LockTime_ );
      #endif

      // If something converged, stick the converged vectors in X
      newlyConverged.resize(newNumConverged);
      for(int i=0; i<newNumConverged; i++)
        newlyConverged[i] = unconvergedIndices[indicesToRemove[i]];

      Teuchos::RCP<MV> helperLocX = MVT::CloneCopy(*locX,indicesToRemove);

      MVT::SetBlock(*helperLocX,newlyConverged,*X_);

      // If everything has converged, we are done
      if(newNumConverged == ncols)
        return;

      // Update unconverged indices
      std::vector<int> helperVec(ncols - newNumConverged);

      std::set_difference(unconvergedIndices.begin(), unconvergedIndices.end(), newlyConverged.begin(), newlyConverged.end(), helperVec.begin());
      unconvergedIndices = helperVec;

      // Get indices of things we want to keep
      std::vector<int> thingsToKeep(ncols - newNumConverged);
      helperVec.resize(ncols);
      for(int i=0; i<ncols; i++)
        helperVec[i] = i;
      ncols = ncols - newNumConverged;

      std::set_difference(helperVec.begin(), helperVec.end(), indicesToRemove.begin(), indicesToRemove.end(), thingsToKeep.begin());

      // Shrink the multivectors
      Teuchos::RCP<MV> helperMV;
      helperMV = MVT::CloneCopy(*V,thingsToKeep);
      V = helperMV;
      helperMV = MVT::CloneCopy(*Y,thingsToKeep);
      Y = helperMV;
      helperMV = MVT::CloneCopy(*R1,thingsToKeep);
      R1 = helperMV;
      helperMV = MVT::CloneCopy(*R2,thingsToKeep);
      R2 = helperMV;
      helperMV = MVT::CloneCopy(*W,thingsToKeep);
      W = helperMV;
      helperMV = MVT::CloneCopy(*W1,thingsToKeep);
      W1 = helperMV;
      helperMV = MVT::CloneCopy(*W2,thingsToKeep);
      W2 = helperMV;
      helperMV = MVT::CloneCopy(*locX,thingsToKeep);
      locX = helperMV;
      scaleHelper = MVT::Clone(*V,ncols);

      // Shrink the arrays
      alpha.resize(ncols);
      oldeps.resize(ncols);
      delta.resize(ncols);
      gbar.resize(ncols);
      phi.resize(ncols);
      gamma.resize(ncols);
      tmpvec.resize(ncols);
      std::vector<Scalar> helperVecS(ncols);
      for(int i=0; i<ncols; i++)
        helperVecS[i] = beta[thingsToKeep[i]];
      beta = helperVecS;
      for(int i=0; i<ncols; i++)
        helperVecS[i] = beta1[thingsToKeep[i]];
      beta1 = helperVecS;
      for(int i=0; i<ncols; i++)
        helperVecS[i] = phibar[thingsToKeep[i]];
      phibar = helperVecS;
      for(int i=0; i<ncols; i++)
        helperVecS[i] = oldBeta[thingsToKeep[i]];
      oldBeta = helperVecS;
      for(int i=0; i<ncols; i++)
        helperVecS[i] = epsln[thingsToKeep[i]];
      epsln = helperVecS;
      for(int i=0; i<ncols; i++)
        helperVecS[i] = cs[thingsToKeep[i]];
      cs = helperVecS;
      for(int i=0; i<ncols; i++)
        helperVecS[i] = sn[thingsToKeep[i]];
      sn = helperVecS;
      for(int i=0; i<ncols; i++)
        helperVecS[i] = dbar[thingsToKeep[i]];
      dbar = helperVecS;

      // Tell operator about the removed indices
      A_->removeIndices(indicesToRemove);
    }

    // Normalize previous vector
    // V = Y / beta
    {
      #ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor lcltimer( *AssignTime_ );
      #endif
      MVT::Assign(*Y,*V);
    }
    for(int i=0; i<ncols; i++)
      tmpvec[i] = ONE/beta[i];
    {
      #ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor lcltimer( *ScaleTime_ );
      #endif
      MVT::MvScale (*V, tmpvec);
    }

    // Apply operator
    // Y = AV
    {
      #ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor lcltimer( *ApplyOpTime_ );
      #endif
      OPT::Apply(*A_, *V, *Y);
    }

    if(iter > 1)
    {
      // Y = Y - beta/oldBeta R1
      for(int i=0; i<ncols; i++)
        tmpvec[i] = beta[i]/oldBeta[i];
      {
        #ifdef ANASAZI_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor lcltimer( *AssignTime_ );
        #endif
        MVT::Assign(*R1,*scaleHelper);
      }
      {
        #ifdef ANASAZI_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor lcltimer( *ScaleTime_ );
        #endif
        MVT::MvScale(*scaleHelper,tmpvec);
      }
      {
        #ifdef ANASAZI_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor lcltimer( *AddTime_ );
        #endif
        MVT::MvAddMv(ONE, *Y, -ONE, *scaleHelper, *Y);
      }
    }

    // alpha = V'*Y
    {
      #ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor lcltimer( *DotTime_ );
      #endif
      MVT::MvDot(*V, *Y, alpha);
    }

    // Y = Y - alpha/beta R2
    for(int i=0; i<ncols; i++)
      tmpvec[i] = alpha[i]/beta[i];
    {
      #ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor lcltimer( *AssignTime_ );
      #endif
      MVT::Assign(*R2,*scaleHelper);
    }    
    {
      #ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor lcltimer( *ScaleTime_ );
      #endif
      MVT::MvScale(*scaleHelper,tmpvec);
    }
    {
      #ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor lcltimer( *AddTime_ );
      #endif
      MVT::MvAddMv(ONE, *Y, -ONE, *scaleHelper, *Y);
    }

    // R1 = R2
    // R2 = Y
    swapHelper = R1;
    R1 = R2;
    R2 = Y;
    Y = swapHelper;

    // Y = M\R2
    if(Prec_ != Teuchos::null)
    {
      #ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor lcltimer( *ApplyPrecTime_ );
      #endif

      OPT::Apply(*Prec_,*R2,*Y);
    }
    else
    {
      #ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor lcltimer( *AssignTime_ );
      #endif
      MVT::Assign(*R2,*Y);
    }

    // Get new beta
    // beta = sqrt(R2'*Y)
    oldBeta = beta;
    {
      #ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor lcltimer( *DotTime_ );
      #endif
      MVT::MvDot(*Y,*R2, beta);

      for(int i=0; i<ncols; i++)
        beta[i] = sqrt(beta[i]);
    }

    // Apply previous rotation
    oldeps = epsln;
    for(int i=0; i<ncols; i++)
    {
      delta[i] = cs[i]*dbar[i]    + sn[i]*alpha[i];
      gbar[i]  = sn[i]*dbar[i]    - cs[i]*alpha[i];
      epsln[i] =                    sn[i]*beta[i];
      dbar[i]  =                  - cs[i]*beta[i];
    }

    // Compute the next plane rotation
    for(int i=0; i<ncols; i++)
    {
      symOrtho(gbar[i], beta[i], &cs[i], &sn[i], &gamma[i]);

      phi[i] = cs[i]*phibar[i];
      phibar[i] = sn[i]*phibar[i];
    }

    // w1 = w2
    // w2 = w
    swapHelper = W1;
    W1 = W2;
    W2 = W;
    W = swapHelper;

    // W = (V - oldeps*W1 - delta*W2) / gamma
    {
      #ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor lcltimer( *AssignTime_ );
      #endif
      MVT::Assign(*W1,*scaleHelper);
    }    
    {
      #ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor lcltimer( *ScaleTime_ );
      #endif
      MVT::MvScale(*scaleHelper,oldeps);
    }
    {
      #ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor lcltimer( *AddTime_ );
      #endif
      MVT::MvAddMv(ONE, *V, -ONE, *scaleHelper, *W);
    }
    {
      #ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor lcltimer( *AssignTime_ );
      #endif
      MVT::Assign(*W2,*scaleHelper);
    }    
    {
      #ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor lcltimer( *ScaleTime_ );
      #endif
      MVT::MvScale(*scaleHelper,delta);
    }
    {
      #ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor lcltimer( *AddTime_ );
      #endif
      MVT::MvAddMv(ONE, *W, -ONE, *scaleHelper, *W);
    }
    for(int i=0; i<ncols; i++)
      tmpvec[i] = ONE/gamma[i];
    {
      #ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor lcltimer( *ScaleTime_ );
      #endif
      MVT::MvScale(*W,tmpvec);
    }

    // Update X
    // X = X + phi*W
    {
      #ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor lcltimer( *AssignTime_ );
      #endif
      MVT::Assign(*W,*scaleHelper);
    }    
    {
      #ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor lcltimer( *ScaleTime_ );
      #endif
      MVT::MvScale(*scaleHelper,phi);
    }
    {
      #ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor lcltimer( *AddTime_ );
      #endif
      MVT::MvAddMv(ONE, *locX, ONE, *scaleHelper, *locX);
    }
  }

  // Stick unconverged vectors in X
  {
    #ifdef ANASAZI_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor lcltimer( *AssignTime_ );
    #endif
    MVT::SetBlock(*locX,unconvergedIndices,*X_);
  }
}

template <class Scalar, class MV, class OP>
void PseudoBlockMinres<Scalar,MV,OP>::symOrtho(Scalar a, Scalar b, Scalar *c, Scalar *s, Scalar *r)
{
  const Scalar absA = std::abs(a);
  const Scalar absB = std::abs(b);
  if ( absB == ZERO ) {
    *s = ZERO;
    *r = absA;
    if ( absA == ZERO )
      *c = ONE;
    else
      *c = a / absA;
  } else if ( absA == ZERO ) {
    *c = ZERO;
    *s = b / absB;
    *r = absB;
  } else if ( absB >= absA ) { // && a!=0 && b!=0
    Scalar tau = a / b;
    if ( b < ZERO )
      *s = -ONE / sqrt( ONE+tau*tau );
    else
      *s =  ONE / sqrt( ONE+tau*tau );
    *c = *s * tau;
    *r = b / *s;
  } else { // (absA > absB) && a!=0 && b!=0
    Scalar tau = b / a;
    if ( a < ZERO )
      *c = -ONE / sqrt( ONE+tau*tau );
    else
      *c =  ONE / sqrt( ONE+tau*tau );
    *s = *c * tau;
    *r = a / *c;
  }
}

}} // End of namespace

#endif
