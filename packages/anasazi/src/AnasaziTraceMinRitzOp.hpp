// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file AnasaziTraceMinRitzOp.hpp
 *  Defines several operators for use with TraceMin
*/

#ifndef TRACEMIN_RITZ_OP_HPP
#define TRACEMIN_RITZ_OP_HPP

#include "AnasaziConfigDefs.hpp"
#include "AnasaziMatOrthoManager.hpp"
#include "AnasaziMinres.hpp"
#include "AnasaziTraceMinBase.hpp"

#ifdef HAVE_ANASAZI_BELOS
  #include "BelosMultiVecTraits.hpp"
  #include "BelosLinearProblem.hpp"
  #include "BelosPseudoBlockGmresSolMgr.hpp"
  #include "BelosOperator.hpp"
  #ifdef HAVE_ANASAZI_TPETRA
    #include "BelosTpetraAdapter.hpp"
  #endif
  #ifdef HAVE_ANASAZI_EPETRA
    #include "BelosEpetraAdapter.hpp"
  #endif
#endif

#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseSolver.hpp"
#include "Teuchos_ScalarTraitsDecl.hpp"


using Teuchos::RCP;

namespace Anasazi {
namespace Experimental {



///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
// Abstract base class for all operators
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
template <class Scalar, class MV, class OP>
class TraceMinOp
{
public:
  virtual ~TraceMinOp() { };
  virtual void Apply(const MV& X, MV& Y) const =0;
  virtual void removeIndices(const std::vector<int>& indicesToRemove) =0;
};



///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
// Defines a projector
// Applies P_i to each individual vector x_i
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
template <class Scalar, class MV, class OP>
class TraceMinProjOp
{
  typedef Anasazi::MultiVecTraits<Scalar,MV>    MVT;
  const Scalar ONE; 

public:
  // Constructors
  TraceMinProjOp(const Teuchos::RCP<const MV> X, const Teuchos::RCP<const OP> B, Teuchos::RCP<Anasazi::MatOrthoManager<Scalar,MV,OP> >  orthman = Teuchos::null);
  TraceMinProjOp(const Teuchos::RCP<const MV> X, const Teuchos::RCP<const OP> B, Teuchos::RCP<Anasazi::MatOrthoManager<Scalar,MV,OP> >  orthman, Teuchos::Array<Teuchos::RCP<const MV> > auxVecs);

  // Destructor
  ~TraceMinProjOp();

  // Applies the projector to a multivector
  void Apply(const MV& X, MV& Y) const;

  // Called by MINRES when certain vectors converge
  void removeIndices(const std::vector<int>& indicesToRemove);

private:  
  Teuchos::Array< Teuchos::RCP<const MV> > projVecs_;
  Teuchos::RCP<Anasazi::MatOrthoManager<Scalar,MV,OP> > orthman_;
  Teuchos::RCP<const OP> B_;

#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
  RCP<Teuchos::Time> ProjTime_;
#endif
};


///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
// This class defines an operator A + \sigma B
// This is used to apply shifts within TraceMin
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
template <class Scalar, class MV, class OP>
class TraceMinRitzOp : public TraceMinOp<Scalar,MV,OP>
{
  template <class Scalar_, class MV_, class OP_> friend class TraceMinProjRitzOp;
  template <class Scalar_, class MV_, class OP_> friend class TraceMinProjRitzOpWithPrec;
  template <class Scalar_, class MV_, class OP_> friend class TraceMinProjectedPrecOp;

  typedef Anasazi::MultiVecTraits<Scalar,MV>     MVT;
  typedef Anasazi::OperatorTraits<Scalar,MV,OP>  OPT;
  const Scalar ONE;  
  const Scalar ZERO;

public:
  // constructors for standard/generalized EVP
  TraceMinRitzOp(const Teuchos::RCP<const OP>& A, const Teuchos::RCP<const OP>& B = Teuchos::null, const Teuchos::RCP<const OP>& Prec = Teuchos::null);

  // Destructor
  ~TraceMinRitzOp() { };

  // sets the Ritz shift
  void setRitzShifts(std::vector<Scalar> shifts) {ritzShifts_ = shifts;};

  Scalar getRitzShift(const int subscript) { return ritzShifts_[subscript]; };

  Teuchos::RCP<TraceMinRitzOp<Scalar,MV,OP> > getPrec() { return Prec_; };

  // sets the tolerances for the inner solves
  void setInnerTol(const std::vector<Scalar>& tolerances) { tolerances_ = tolerances; };

  void getInnerTol(std::vector<Scalar>& tolerances) const { tolerances = tolerances_; };

  void setMaxIts(const int maxits) { maxits_ = maxits; };
  
  int getMaxIts() const { return maxits_; };

  // applies A+\sigma B to a vector
  void Apply(const MV& X, MV& Y) const;

  // returns (A+\sigma B)\X
  void ApplyInverse(const MV& X, MV& Y);

  void removeIndices(const std::vector<int>& indicesToRemove);

private:  
  Teuchos::RCP<const OP> A_, B_;
  Teuchos::RCP<TraceMinRitzOp<Scalar,MV,OP> > Prec_;

  int maxits_;
  std::vector<Scalar> ritzShifts_;
  std::vector<Scalar> tolerances_;

  Teuchos::RCP< PseudoBlockMinres< Scalar,MV,TraceMinRitzOp<Scalar,MV,OP> > > solver_;

#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
  RCP<Teuchos::Time> PetraMultTime_, AopMultTime_;
#endif
};



///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
// Defines an operator P (A + \sigma B) P
// Used for TraceMin with the projected iterative solver
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
template <class Scalar, class MV, class OP>
class TraceMinProjRitzOp : public TraceMinOp<Scalar,MV,OP>
{
  typedef Anasazi::MultiVecTraits<Scalar,MV>                                MVT;
  typedef Anasazi::OperatorTraits<Scalar,MV,TraceMinRitzOp<Scalar,MV,OP> >  OPT;

public:
  // constructors for standard/generalized EVP
  TraceMinProjRitzOp(const Teuchos::RCP<TraceMinRitzOp<Scalar,MV,OP> >& Op, const Teuchos::RCP<const MV> projVecs, const Teuchos::RCP<Anasazi::MatOrthoManager<Scalar,MV,OP> > orthman = Teuchos::null);
  TraceMinProjRitzOp(const Teuchos::RCP<TraceMinRitzOp<Scalar,MV,OP> >& Op, const Teuchos::RCP<const MV> projVecs, const Teuchos::RCP<Anasazi::MatOrthoManager<Scalar,MV,OP> > orthman, Teuchos::Array<Teuchos::RCP<const MV> > auxVecs);

  // applies P (A+\sigma B) P to a vector
  void Apply(const MV& X, MV& Y) const;

  // returns P(A+\sigma B)P\X
  // this is not const due to the clumsiness with amSolving
  void ApplyInverse(const MV& X, MV& Y);

  void removeIndices(const std::vector<int>& indicesToRemove);

private:  
  Teuchos::RCP< TraceMinRitzOp<Scalar,MV,OP> > Op_;
  Teuchos::RCP< TraceMinProjOp<Scalar,MV,OP> > projector_;

  Teuchos::RCP< PseudoBlockMinres< Scalar,MV,TraceMinProjRitzOp<Scalar,MV,OP> > > solver_;
};



///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
// Defines a preconditioner to be used with our projection method
// Because we're using projected CG/minres/gmres, this preconditioner has to do projection as well
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
// TODO: Should this be public?
template <class Scalar, class MV, class OP>
class TraceMinProjectedPrecOp : public TraceMinOp<Scalar,MV,OP>
{
  typedef Anasazi::MultiVecTraits<Scalar,MV>     MVT;
  typedef Anasazi::OperatorTraits<Scalar,MV,OP>  OPT;
  const Scalar ONE; 

public:
  // constructors for standard/generalized EVP
  TraceMinProjectedPrecOp(const Teuchos::RCP<const OP> Op, const Teuchos::RCP<const MV> projVecs, Teuchos::RCP< Anasazi::MatOrthoManager<Scalar,MV,OP> > orthman = Teuchos::null);
  TraceMinProjectedPrecOp(const Teuchos::RCP<const OP> Op, const Teuchos::RCP<const MV> projVecs, Teuchos::RCP< Anasazi::MatOrthoManager<Scalar,MV,OP> > orthman, Teuchos::Array< Teuchos::RCP<const MV> > auxVecs);

  ~TraceMinProjectedPrecOp();

  void Apply(const MV& X, MV& Y) const;

  void removeIndices(const std::vector<int>& indicesToRemove);

private:  
  Teuchos::RCP<const OP> Op_;
  Teuchos::Array< Teuchos::RCP<const MV> > projVecs_;

  Teuchos::RCP< Anasazi::MatOrthoManager<Scalar,MV,OP> > orthman_;
  Teuchos::RCP<const OP> B_;
};



///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
// Defines a preconditioner to be used with our projection method
// Because we're using projected CG/minres/gmres, this preconditioner has to do projection as well
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef HAVE_ANASAZI_BELOS
template <class Scalar, class MV, class OP>
class TraceMinProjRitzOpWithPrec : public TraceMinOp<Scalar,MV,OP>
{
  typedef Anasazi::MultiVecTraits<Scalar,MV>                                MVT;
  typedef Anasazi::OperatorTraits<Scalar,MV,TraceMinRitzOp<Scalar,MV,OP> >  OPT;
  const Scalar ONE; 

public:
  // constructors for standard/generalized EVP
  TraceMinProjRitzOpWithPrec(const Teuchos::RCP<TraceMinRitzOp<Scalar,MV,OP> >& Op, const Teuchos::RCP<const MV> projVecs, Teuchos::RCP< Anasazi::MatOrthoManager<Scalar,MV,OP> > orthman = Teuchos::null);
  TraceMinProjRitzOpWithPrec(const Teuchos::RCP<TraceMinRitzOp<Scalar,MV,OP> >& Op, const Teuchos::RCP<const MV> projVecs, Teuchos::RCP< Anasazi::MatOrthoManager<Scalar,MV,OP> > orthman, Teuchos::Array< Teuchos::RCP<const MV> > auxVecs);

  void Apply(const MV& X, MV& Y) const;

  // returns P(A+\sigma B)P\X
  // this is not const due to the clumsiness with amSolving
  void ApplyInverse(const MV& X, MV& Y);

  void removeIndices(const std::vector<int>& indicesToRemove);

private:  
  Teuchos::RCP<TraceMinRitzOp<Scalar,MV,OP> > Op_;
  Teuchos::RCP<TraceMinProjectedPrecOp<Scalar,MV,OP> > Prec_;

  Teuchos::RCP< Belos::PseudoBlockGmresSolMgr< Scalar,MV,TraceMinOp<Scalar,MV,OP> > > solver_;
  Teuchos::RCP< Belos::LinearProblem<Scalar,MV,TraceMinOp<Scalar,MV,OP> > > problem_;
};
#endif

}} // end of namespace



///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
// Operator traits classes
// Required to use user-defined operators with a Krylov solver in Belos
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
namespace Anasazi
{
template <class Scalar, class MV, class OP>
class OperatorTraits <Scalar, MV, Experimental::TraceMinOp<Scalar,MV,OP> >
  {
  public:
    static void
    Apply (const Experimental::TraceMinOp<Scalar,MV,OP>& Op,
           const MV& x,
           MV& y) {Op.Apply(x,y);};

    //! Whether Op implements applying the transpose.
    static bool
    HasApplyTranspose (const Experimental::TraceMinOp<Scalar,MV,OP>& Op) {return true;};
  };
}


#ifdef HAVE_ANASAZI_BELOS
namespace Belos
{
template <class Scalar, class MV, class OP>
class OperatorTraits <Scalar, MV, Anasazi::Experimental::TraceMinOp<Scalar,MV,OP> >
  {
  public:
    static void
    Apply (const Anasazi::Experimental::TraceMinOp<Scalar,MV,OP>& Op,
           const MV& x,
           MV& y,  Belos::ETrans trans = NOTRANS) {Op.Apply(x,y);};

    //! Whether Op implements applying the transpose.
    static bool
    HasApplyTranspose (const Anasazi::Experimental::TraceMinOp<Scalar,MV,OP>& Op) {return true;};
  };
}
#endif



namespace Anasazi
{
template <class Scalar, class MV, class OP>
class OperatorTraits <Scalar, MV, Experimental::TraceMinRitzOp<Scalar,MV,OP> >
  {
  public:
    static void
    Apply (const Experimental::TraceMinRitzOp<Scalar,MV,OP>& Op,
           const MV& x,
           MV& y) {Op.Apply(x,y);};

    //! Whether Op implements applying the transpose.
    static bool
    HasApplyTranspose (const Experimental::TraceMinRitzOp<Scalar,MV,OP>& Op) {return true;};
  };
}



namespace Anasazi
{
template <class Scalar, class MV, class OP>
class OperatorTraits <Scalar, MV, Experimental::TraceMinProjRitzOp<Scalar,MV,OP> >
  {
  public:
    static void
    Apply (const Experimental::TraceMinProjRitzOp<Scalar,MV,OP>& Op,
           const MV& x,
           MV& y) {Op.Apply(x,y);};

    //! Whether Op implements applying the transpose.
    static bool
    HasApplyTranspose (const Experimental::TraceMinProjRitzOp<Scalar,MV,OP>& Op) {return true;};
  };
}



namespace Anasazi
{
template <class Scalar, class MV, class OP>
class OperatorTraits <Scalar, MV, Experimental::TraceMinProjectedPrecOp<Scalar,MV,OP> >
  {
  public:
    static void
    Apply (const Experimental::TraceMinProjectedPrecOp<Scalar,MV,OP>& Op,
           const MV& x,
           MV& y) {Op.Apply(x,y);};

    //! Whether Op implements applying the transpose.
    static bool
    HasApplyTranspose (const Experimental::TraceMinProjectedPrecOp<Scalar,MV,OP>& Op) {return true;};
  };
}



namespace Anasazi {
namespace Experimental {
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
// TraceMinProjOp implementations
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
template <class Scalar, class MV, class OP>
TraceMinProjOp<Scalar,MV,OP>::TraceMinProjOp(const Teuchos::RCP<const MV> X, const Teuchos::RCP<const OP> B, Teuchos::RCP<Anasazi::MatOrthoManager<Scalar,MV,OP> > orthman) :
ONE(Teuchos::ScalarTraits<Scalar>::one())
{
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
  ProjTime_ = Teuchos::TimeMonitor::getNewTimer("Anasazi: TraceMinProjOp::Apply()");
#endif

  B_ = B;
  orthman_ = orthman;
  if(orthman_ != Teuchos::null && B_ != Teuchos::null)
    orthman_->setOp(Teuchos::null);

  // Make it so X'BBX = I
  // If there is no B, this step is unnecessary because X'X = I already
  if(B_ != Teuchos::null)
  {
    int nvec = MVT::GetNumberVecs(*X);

    if(orthman_ != Teuchos::null)
    {
      // Want: X'X = I NOT X'MX = I
      Teuchos::RCP<MV> helperMV = MVT::CloneCopy(*X);
      orthman_->normalizeMat(*helperMV);
      projVecs_.push_back(helperMV);
    }
    else
    {
      std::vector<Scalar> normvec(nvec);
      MVT::MvNorm(*X,normvec);
      for(int i=0; i<nvec; i++)
        normvec[i] = ONE/normvec[i];
      Teuchos::RCP<MV> helperMV = MVT::CloneCopy(*X);
      MVT::MvScale(*helperMV,normvec);
      projVecs_.push_back(helperMV);
    }
  }
  else
    projVecs_.push_back(MVT::CloneCopy(*X));
}


template <class Scalar, class MV, class OP>
TraceMinProjOp<Scalar,MV,OP>::TraceMinProjOp(const Teuchos::RCP<const MV> X, const Teuchos::RCP<const OP> B, Teuchos::RCP<Anasazi::MatOrthoManager<Scalar,MV,OP> >  orthman, Teuchos::Array<Teuchos::RCP<const MV> > auxVecs) :
ONE(Teuchos::ScalarTraits<Scalar>::one())
{
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
  ProjTime_ = Teuchos::TimeMonitor::getNewTimer("Anasazi: TraceMinProjOp::Apply()");
#endif

  B_ = B;
  orthman_ = orthman;
  if(B_ != Teuchos::null)
    orthman_->setOp(Teuchos::null);

  projVecs_ = auxVecs;

  // Make it so X'BBX = I
  // If there is no B, this step is unnecessary because X'X = I already
  if(B_ != Teuchos::null)
  {
    // Want: X'X = I NOT X'MX = I
    Teuchos::RCP<MV> helperMV = MVT::CloneCopy(*X);
    orthman_->normalizeMat(*helperMV);
    projVecs_.push_back(helperMV);

  }
  else
    projVecs_.push_back(MVT::CloneCopy(*X));
}


// Destructor - make sure to reset the operator in the ortho manager
template <class Scalar, class MV, class OP>
TraceMinProjOp<Scalar,MV,OP>::~TraceMinProjOp()
{
  if(orthman_ != Teuchos::null)
    orthman_->setOp(B_);
}


// Compute Px = x - proj proj'x
template <class Scalar, class MV, class OP>
void TraceMinProjOp<Scalar,MV,OP>::Apply(const MV& X, MV& Y) const
{
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor lcltimer( *ProjTime_ );
#endif

  if(orthman_ != Teuchos::null)
  {
    MVT::Assign(X,Y);
    orthman_->projectMat(Y,projVecs_);
  }
  else
  {
    int nvec = MVT::GetNumberVecs(X);
    std::vector<Scalar> dotProducts(nvec);
    MVT::MvDot(*projVecs_[0],X,dotProducts);
    Teuchos::RCP<MV> helper = MVT::CloneCopy(*projVecs_[0]);
    MVT::MvScale(*helper,dotProducts);
    MVT::MvAddMv(ONE,X,-ONE,*helper,Y);
  }
}



template <class Scalar, class MV, class OP>
void TraceMinProjOp<Scalar,MV,OP>::removeIndices(const std::vector<int>& indicesToRemove)
{
  if (orthman_ == Teuchos::null) {
    const int nprojvecs = projVecs_.size();
    const int nvecs = MVT::GetNumberVecs(*projVecs_[nprojvecs-1]);
    const int numRemoving = indicesToRemove.size();
    std::vector<int> helper(nvecs), indicesToLeave(nvecs-numRemoving);

    for (int i=0; i<nvecs; i++) {
      helper[i] = i;
    }

    std::set_difference(helper.begin(), helper.end(), indicesToRemove.begin(), indicesToRemove.end(), indicesToLeave.begin());

    Teuchos::RCP<MV> helperMV = MVT::CloneCopy(*projVecs_[nprojvecs-1],indicesToLeave);
    projVecs_[nprojvecs-1] = helperMV;
  }
}


///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
// TraceMinRitzOp implementations
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
template <class Scalar, class MV, class OP>
TraceMinRitzOp<Scalar,MV,OP>::TraceMinRitzOp(const Teuchos::RCP<const OP>& A, const Teuchos::RCP<const OP>& B, const Teuchos::RCP<const OP>& Prec) :
ONE(Teuchos::ScalarTraits<Scalar>::one()),
ZERO(Teuchos::ScalarTraits<Scalar>::zero())
{
  A_ = A;
  B_ = B;
  // TODO: maxits should not be hard coded
  maxits_ = 200;

#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
  PetraMultTime_ = Teuchos::TimeMonitor::getNewTimer("Anasazi: TraceMinRitzOp: *Petra::Apply()");
  AopMultTime_ = Teuchos::TimeMonitor::getNewTimer("Anasazi: TraceMinRitzOp::Apply()");
#endif

  // create the operator for my minres solver
  Teuchos::RCP<TraceMinRitzOp<Scalar,MV,OP> > linSolOp = Teuchos::rcp(this);
  linSolOp.release();

  // TODO: This should support left and right prec
  if(Prec != Teuchos::null)    
    Prec_ = Teuchos::rcp( new TraceMinRitzOp<Scalar,MV,OP>(Prec) );
  
  // create my minres solver
  solver_ = Teuchos::rcp( new PseudoBlockMinres< Scalar,MV,TraceMinRitzOp<Scalar,MV,OP> >(linSolOp,Prec_) );
}



template <class Scalar, class MV, class OP>
void TraceMinRitzOp<Scalar,MV,OP>::Apply(const MV& X, MV& Y) const
{
  int nvecs = MVT::GetNumberVecs(X);

#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor outertimer( *AopMultTime_ );
#endif

  // Y := A*X
  {
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor lcltimer( *PetraMultTime_ );
#endif

    OPT::Apply(*A_,X,Y);
  }

  // If we are a preconditioner, we're not using shifts
  if(ritzShifts_.size() > 0)
  {
    // Get the indices of nonzero Ritz shifts
    std::vector<int> nonzeroRitzIndices;
    nonzeroRitzIndices.reserve(nvecs);
    for(int i=0; i<nvecs; i++)
    {
      if(ritzShifts_[i] != ZERO)
        nonzeroRitzIndices.push_back(i);
    }

    // Handle Ritz shifts
    int numRitzShifts = nonzeroRitzIndices.size();
    if(numRitzShifts > 0)
    {
      // Get pointers to the appropriate parts of X and Y
      Teuchos::RCP<const MV> ritzX = MVT::CloneView(X,nonzeroRitzIndices);
      Teuchos::RCP<MV> ritzY = MVT::CloneViewNonConst(Y,nonzeroRitzIndices);

      // Get the nonzero Ritz shifts
      std::vector<Scalar> nonzeroRitzShifts(numRitzShifts);
      for(int i=0; i<numRitzShifts; i++)
        nonzeroRitzShifts[i] = ritzShifts_[nonzeroRitzIndices[i]];

      // Compute Y := AX + ritzShift BX
      if(B_ != Teuchos::null)
      {
        Teuchos::RCP<MV> BX = MVT::Clone(*ritzX,numRitzShifts);
        OPT::Apply(*B_,*ritzX,*BX);
        MVT::MvScale(*BX,nonzeroRitzShifts);
        MVT::MvAddMv(ONE,*ritzY,-ONE,*BX,*ritzY);
      }
      // Compute Y := AX + ritzShift X
      else
      {
        Teuchos::RCP<MV> scaledX = MVT::CloneCopy(*ritzX);
        MVT::MvScale(*scaledX,nonzeroRitzShifts);
        MVT::MvAddMv(ONE,*ritzY,-ONE,*scaledX,*ritzY);
      }
    }
  }
}



template <class Scalar, class MV, class OP>
void TraceMinRitzOp<Scalar,MV,OP>::ApplyInverse(const MV& X, MV& Y)
{
  int nvecs = MVT::GetNumberVecs(X);
  std::vector<int> indices(nvecs);
  for(int i=0; i<nvecs; i++)
    indices[i] = i;

  Teuchos::RCP<const MV> rcpX = MVT::CloneView(X,indices);
  Teuchos::RCP<MV> rcpY = MVT::CloneViewNonConst(Y,indices);

  // Solve the linear system A*Y = X
  solver_->setTol(tolerances_);
  solver_->setMaxIter(maxits_);

  // Set solution and RHS
  solver_->setSol(rcpY);
  solver_->setRHS(rcpX);

  // Solve the linear system
  solver_->solve();  
}



template <class Scalar, class MV, class OP>
void TraceMinRitzOp<Scalar,MV,OP>::removeIndices(const std::vector<int>& indicesToRemove)
{
  int nvecs = tolerances_.size();
  int numRemoving = indicesToRemove.size();
  std::vector<int> helper(nvecs), indicesToLeave(nvecs-numRemoving);
  std::vector<Scalar> helperS(nvecs-numRemoving);

  for(int i=0; i<nvecs; i++)
    helper[i] = i;

  std::set_difference(helper.begin(), helper.end(), indicesToRemove.begin(), indicesToRemove.end(), indicesToLeave.begin());

  for(int i=0; i<nvecs-numRemoving; i++)
    helperS[i] = ritzShifts_[indicesToLeave[i]];
  ritzShifts_ = helperS;

  for(int i=0; i<nvecs-numRemoving; i++)
    helperS[i] = tolerances_[indicesToLeave[i]];
  tolerances_ = helperS;
}



///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
// TraceMinProjRitzOp implementations
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
template <class Scalar, class MV, class OP>
TraceMinProjRitzOp<Scalar,MV,OP>::TraceMinProjRitzOp(const Teuchos::RCP<TraceMinRitzOp<Scalar,MV,OP> >& Op, const Teuchos::RCP<const MV> projVecs, const Teuchos::RCP<Anasazi::MatOrthoManager<Scalar,MV,OP> > orthman) 
{ 
  Op_ = Op; 

  // Create the projector object
  projector_ = Teuchos::rcp( new TraceMinProjOp<Scalar,MV,OP>(projVecs, Op_->B_, orthman) );

  // create the operator for my minres solver
  Teuchos::RCP<TraceMinProjRitzOp<Scalar,MV,OP> > linSolOp = Teuchos::rcp(this);
  linSolOp.release();

  // create my minres solver
  solver_ = Teuchos::rcp( new PseudoBlockMinres< Scalar,MV,TraceMinProjRitzOp<Scalar,MV,OP> >(linSolOp) );
}


template <class Scalar, class MV, class OP>
TraceMinProjRitzOp<Scalar,MV,OP>::TraceMinProjRitzOp(const Teuchos::RCP<TraceMinRitzOp<Scalar,MV,OP> >& Op, const Teuchos::RCP<const MV> projVecs, const Teuchos::RCP<Anasazi::MatOrthoManager<Scalar,MV,OP> > orthman, Teuchos::Array<Teuchos::RCP<const MV> > auxVecs)
{
  Op_ = Op; 

  // Create the projector object
  projector_ = Teuchos::rcp( new TraceMinProjOp<Scalar,MV,OP>(projVecs, Op_->B_, orthman, auxVecs) );

  // create the operator for my minres solver
  Teuchos::RCP<TraceMinProjRitzOp<Scalar,MV,OP> > linSolOp = Teuchos::rcp(this);
  linSolOp.release();

  // create my minres solver
  solver_ = Teuchos::rcp( new PseudoBlockMinres< Scalar,MV,TraceMinProjRitzOp<Scalar,MV,OP> >(linSolOp) );
}



// Y:= P (A+\sigma B) P X
template <class Scalar, class MV, class OP>
void TraceMinProjRitzOp<Scalar,MV,OP>::Apply(const MV& X, MV& Y) const
{
  int nvecs = MVT::GetNumberVecs(X);

//  // compute PX
//  Teuchos::RCP<MV> PX = MVT::Clone(X,nvecs);
//  projector_->Apply(X,*PX);

  // compute (A+\sigma B) P X
  Teuchos::RCP<MV> APX = MVT::Clone(X,nvecs);
  OPT::Apply(*Op_,X,*APX);

  // compute Y := P (A+\sigma B) P X
  projector_->Apply(*APX,Y);
}



template <class Scalar, class MV, class OP>
void TraceMinProjRitzOp<Scalar,MV,OP>::ApplyInverse(const MV& X, MV& Y)
{
  int nvecs = MVT::GetNumberVecs(X);
  std::vector<int> indices(nvecs);
  for(int i=0; i<nvecs; i++)
    indices[i] = i;

  Teuchos::RCP<MV> rcpY = MVT::CloneViewNonConst(Y,indices);
  Teuchos::RCP<MV> PX = MVT::Clone(X,nvecs);
  projector_->Apply(X,*PX);

  // Solve the linear system A*Y = X
  solver_->setTol(Op_->tolerances_);
  solver_->setMaxIter(Op_->maxits_);

  // Set solution and RHS
  solver_->setSol(rcpY);
  solver_->setRHS(PX);

  // Solve the linear system
  solver_->solve();  
}



template <class Scalar, class MV, class OP>
void TraceMinProjRitzOp<Scalar,MV,OP>::removeIndices(const std::vector<int>& indicesToRemove)
{
  Op_->removeIndices(indicesToRemove);

  projector_->removeIndices(indicesToRemove);
}




#ifdef HAVE_ANASAZI_BELOS
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
// TraceMinProjRitzOpWithPrec implementations
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
template <class Scalar, class MV, class OP>
TraceMinProjRitzOpWithPrec<Scalar,MV,OP>::TraceMinProjRitzOpWithPrec(const Teuchos::RCP<TraceMinRitzOp<Scalar,MV,OP> >& Op, const Teuchos::RCP<const MV> projVecs, Teuchos::RCP< Anasazi::MatOrthoManager<Scalar,MV,OP> > orthman) :
ONE(Teuchos::ScalarTraits<Scalar>::one())
{
  Op_ = Op;

  // create the operator for the Belos solver
  Teuchos::RCP<TraceMinProjRitzOpWithPrec<Scalar,MV,OP> > linSolOp = Teuchos::rcp(this);
  linSolOp.release();

  // Create the linear problem
  problem_ = Teuchos::rcp(new Belos::LinearProblem< Scalar,MV,TraceMinOp<Scalar,MV,OP> >());

  // Set the operator
  problem_->setOperator(linSolOp);

  // Set the preconditioner
  // TODO: This does not support right preconditioning
  Prec_ = Teuchos::rcp( new TraceMinProjectedPrecOp<Scalar,MV,OP>(Op_->Prec_->A_, projVecs, orthman) );
//  problem_->setLeftPrec(Prec_);

  // create the pseudoblock gmres solver
  // minres has trouble with the projected preconditioner
  solver_ = Teuchos::rcp( new Belos::PseudoBlockGmresSolMgr< Scalar,MV,TraceMinOp<Scalar,MV,OP> >());
}


template <class Scalar, class MV, class OP>
TraceMinProjRitzOpWithPrec<Scalar,MV,OP>::TraceMinProjRitzOpWithPrec(const Teuchos::RCP<TraceMinRitzOp<Scalar,MV,OP> >& Op, const Teuchos::RCP<const MV> projVecs, Teuchos::RCP< Anasazi::MatOrthoManager<Scalar,MV,OP> > orthman, Teuchos::Array< Teuchos::RCP<const MV> > auxVecs) :
ONE(Teuchos::ScalarTraits<Scalar>::one())
{
  Op_ = Op;

  // create the operator for the Belos solver
  Teuchos::RCP<const TraceMinProjRitzOpWithPrec<Scalar,MV,OP> > linSolOp = Teuchos::rcp(this);
  linSolOp.release();

  // Create the linear problem
  problem_ = Teuchos::rcp(new Belos::LinearProblem< Scalar,MV,TraceMinOp<Scalar,MV,OP> >());

  // Set the operator
  problem_->setOperator(linSolOp);

  // Set the preconditioner
  // TODO: This does not support right preconditioning
  Prec_ = Teuchos::rcp( new TraceMinProjectedPrecOp<Scalar,MV,OP>(Op_->Prec_->A_, projVecs, orthman, auxVecs) );
//  problem_->setLeftPrec(Prec_);

  // create the pseudoblock gmres solver
  // minres has trouble with the projected preconditioner
  solver_ = Teuchos::rcp( new Belos::PseudoBlockGmresSolMgr< Scalar,MV,TraceMinOp<Scalar,MV,OP> >());
}


template <class Scalar, class MV, class OP>
void TraceMinProjRitzOpWithPrec<Scalar,MV,OP>::Apply(const MV& X, MV& Y) const
{
  int nvecs = MVT::GetNumberVecs(X);
  RCP<MV> Ydot = MVT::Clone(Y,nvecs);
  OPT::Apply(*Op_,X,*Ydot);
  Prec_->Apply(*Ydot,Y);
}


template <class Scalar, class MV, class OP>
void TraceMinProjRitzOpWithPrec<Scalar,MV,OP>::ApplyInverse(const MV& X, MV& Y)
{
  int nvecs = MVT::GetNumberVecs(X);
  std::vector<int> indices(nvecs);
  for(int i=0; i<nvecs; i++)
    indices[i] = i;

  Teuchos::RCP<MV> rcpY = MVT::CloneViewNonConst(Y,indices);
  Teuchos::RCP<MV> rcpX = MVT::Clone(X,nvecs);
  
  Prec_->Apply(X,*rcpX);
  
  // Create the linear problem
  problem_->setProblem(rcpY,rcpX);

  // Set the problem for the solver
  solver_->setProblem(problem_);

  // Set up the parameters for gmres
  // TODO: Accept maximum number of iterations
  // TODO: Make sure single shift really means single shift
  // TODO: Look into fixing my problem with the deflation quorum
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::rcp(new Teuchos::ParameterList());
  pl->set("Convergence Tolerance", Op_->tolerances_[0]);
  pl->set("Block Size", nvecs);
//  pl->set("Verbosity", Belos::IterationDetails + Belos::StatusTestDetails + Belos::Debug);
//  pl->set("Output Frequency", 1);
  pl->set("Maximum Iterations", Op_->getMaxIts());
  pl->set("Num Blocks", Op_->getMaxIts());
  solver_->setParameters(pl);

  // Solve the linear system
  solver_->solve();  
}


template <class Scalar, class MV, class OP>
void TraceMinProjRitzOpWithPrec<Scalar,MV,OP>::removeIndices(const std::vector<int>& indicesToRemove)
{
  Op_->removeIndices(indicesToRemove);

  Prec_->removeIndices(indicesToRemove);
}
#endif


///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
// TraceMinProjectedPrecOp implementations
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
template <class Scalar, class MV, class OP>
TraceMinProjectedPrecOp<Scalar,MV,OP>::TraceMinProjectedPrecOp(const Teuchos::RCP<const OP> Op, const Teuchos::RCP<const MV> projVecs, Teuchos::RCP< Anasazi::MatOrthoManager<Scalar,MV,OP> > orthman) :
ONE (Teuchos::ScalarTraits<Scalar>::one())
{
  Op_ = Op;
  orthman_ = orthman;

  int nvecs = MVT::GetNumberVecs(*projVecs);
  Teuchos::RCP<MV> helperMV = MVT::Clone(*projVecs,nvecs);

  // Compute Prec \ projVecs
  OPT::Apply(*Op_,*projVecs,*helperMV);

  if(orthman_ != Teuchos::null)
  {
    // Set the operator for the inner products
    B_ = orthman_->getOp();
    orthman_->setOp(Op_);

    Teuchos::RCP<MV> locProjVecs = MVT::CloneCopy(*projVecs);

    // Normalize the vectors such that Y' Prec \ Y = I
    const int rank = orthman_->normalizeMat (*locProjVecs, Teuchos::null, helperMV);

    // FIXME (mfh 08 Aug 2014) Write a named exception for this case.
    TEUCHOS_TEST_FOR_EXCEPTION(rank != nvecs, std::runtime_error, "Belos::TraceMinProjectedPrecOp: constructor: Loss of orthogonality detected.");

    orthman_->setOp(Teuchos::null);
  }
  else
  {
    std::vector<Scalar> dotprods(nvecs);
    MVT::MvDot(*projVecs,*helperMV,dotprods);

    for(int i=0; i<nvecs; i++)
      dotprods[i] = ONE/sqrt(dotprods[i]);

    MVT::MvScale(*helperMV, dotprods);
  }
  
  projVecs_.push_back(helperMV);
}


template <class Scalar, class MV, class OP>
TraceMinProjectedPrecOp<Scalar,MV,OP>::TraceMinProjectedPrecOp(const Teuchos::RCP<const OP> Op, const Teuchos::RCP<const MV> projVecs, Teuchos::RCP< Anasazi::MatOrthoManager<Scalar,MV,OP> > orthman, Teuchos::Array< Teuchos::RCP<const MV> > auxVecs) :
ONE(Teuchos::ScalarTraits<Scalar>::one())
{
  Op_ = Op;
  orthman_ = orthman;

  int nvecs;
  Teuchos::RCP<MV> locProjVecs;

  // Add the aux vecs to the projector
  if(auxVecs.size() > 0)
  {
    // Get the total number of vectors
    nvecs = MVT::GetNumberVecs(*projVecs);
    for(int i=0; i<auxVecs.size(); i++)
      nvecs += MVT::GetNumberVecs(*auxVecs[i]);

    // Allocate space for all of them
    locProjVecs = MVT::Clone(*projVecs, nvecs);

    // Copy the vectors over
    int startIndex = 0;
    std::vector<int> locind(nvecs);

    locind.resize(MVT::GetNumberVecs(*projVecs));
    for (size_t i = 0; i<locind.size(); i++) {
      locind[i] = startIndex + i;
    }
    startIndex += locind.size();
    MVT::SetBlock(*projVecs,locind,*locProjVecs);

    for (size_t i=0; i < static_cast<size_t> (auxVecs.size ()); ++i)
    {
      locind.resize(MVT::GetNumberVecs(*auxVecs[i]));
      for(size_t j=0; j<locind.size(); j++) locind[j] = startIndex + j;
      startIndex += locind.size();
      MVT::SetBlock(*auxVecs[i],locind,*locProjVecs);
    }
  }
  else
  {
    // Copy the vectors over
    nvecs = MVT::GetNumberVecs(*projVecs);
    locProjVecs = MVT::CloneCopy(*projVecs);
  }
  
  Teuchos::RCP<MV> helperMV = MVT::Clone(*projVecs,nvecs);

  // Compute Prec \ projVecs
  OPT::Apply(*Op_,*locProjVecs,*helperMV);
  
  // Set the operator for the inner products
  B_ = orthman_->getOp();
  orthman_->setOp(Op_);

  // Normalize the vectors such that Y' Prec \ Y = I
  const int rank = orthman_->normalizeMat(*locProjVecs,Teuchos::null,helperMV);
  
  projVecs_.push_back(helperMV);

//  helperMV->describe(*(Teuchos::VerboseObjectBase::getDefaultOStream()),Teuchos::VERB_EXTREME);

  // FIXME (mfh 08 Aug 2014) Write a named exception for this case.  
  TEUCHOS_TEST_FOR_EXCEPTION(
          rank != nvecs, std::runtime_error, 
    "Belos::TraceMinProjectedPrecOp: constructor: Loss of orthogonality detected.");

  orthman_->setOp(Teuchos::null);
}


template <class Scalar, class MV, class OP>
TraceMinProjectedPrecOp<Scalar,MV,OP>::~TraceMinProjectedPrecOp()
{
  if(orthman_ != Teuchos::null)
    orthman_->setOp(B_);
}



template <class Scalar, class MV, class OP>
void TraceMinProjectedPrecOp<Scalar,MV,OP>::Apply(const MV& X, MV& Y) const
{
  int nvecsX = MVT::GetNumberVecs(X);

  if(orthman_ != Teuchos::null)
  {
    // Y = M\X - proj proj' X
    int nvecsP = MVT::GetNumberVecs(*projVecs_[0]);
    OPT::Apply(*Op_,X,Y);
	
    Teuchos::RCP< Teuchos::SerialDenseMatrix<int,Scalar> > projX = Teuchos::rcp(new Teuchos::SerialDenseMatrix<int,Scalar>(nvecsP,nvecsX));

    MVT::MvTransMv(ONE, *projVecs_[0], X, *projX);
	
    MVT::MvTimesMatAddMv(-ONE, *projVecs_[0], *projX, ONE, Y);
  }
  else
  {
    Teuchos::RCP<MV> MX = MVT::Clone(X,nvecsX);
    OPT::Apply(*Op_,X,*MX);
	
    std::vector<Scalar> dotprods(nvecsX);
    MVT::MvDot(*projVecs_[0], X, dotprods);
	
    Teuchos::RCP<MV> helper = MVT::CloneCopy(*projVecs_[0]);
    MVT::MvScale(*helper, dotprods);
    MVT::MvAddMv(ONE, *MX, -ONE, *helper, Y);
  }
}

  
template <class Scalar, class MV, class OP>
void TraceMinProjectedPrecOp<Scalar,MV,OP>::removeIndices(const std::vector<int>& indicesToRemove)
{
  if(orthman_ == Teuchos::null)
  {
    int nvecs = MVT::GetNumberVecs(*projVecs_[0]);
    int numRemoving = indicesToRemove.size();
    std::vector<int> helper(nvecs), indicesToLeave(nvecs-numRemoving);

    for(int i=0; i<nvecs; i++)
      helper[i] = i;

    std::set_difference(helper.begin(), helper.end(), indicesToRemove.begin(), indicesToRemove.end(), indicesToLeave.begin());

    Teuchos::RCP<const MV> helperMV = MVT::CloneCopy(*projVecs_[0],indicesToLeave);
    projVecs_[0] = helperMV;
  }
}

}} // end of namespace

#endif
