// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// @HEADER

#ifndef ANASAZI_STATUS_TEST_SPECTRANS_HPP
#define ANASAZI_STATUS_TEST_SPECTRANS_HPP

#include "AnasaziTypes.hpp"
#include "AnasaziStatusTestResNorm.hpp"
#include "AnasaziTraceMinBase.hpp"

using Teuchos::RCP;

namespace Anasazi {
namespace Experimental {

  template<class ScalarType, class MV, class OP>
  class StatusTestSpecTrans : public StatusTest<ScalarType,MV,OP> {

  typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
  typedef MultiVecTraits<ScalarType,MV>                             MVT;
  typedef OperatorTraits<ScalarType,MV,OP>                          OPT;

  public:

    // Constructor
    StatusTestSpecTrans(MagnitudeType tol, int quorum = -1, ResType whichNorm = RES_2NORM, bool scaled = true, bool throwExceptionOnNan = true, const RCP<const OP> Mop = Teuchos::null);

    // Destructor
    virtual ~StatusTestSpecTrans() {};

    // Check whether the test passed or failed
    TestStatus checkStatus(Eigensolver<ScalarType,MV,OP> *solver);

    // Return the result of the most recent checkStatus call
    TestStatus getStatus() const { return state_; }

    // Get the indices for the vectors that passed the test
    std::vector<int> whichVecs() const { return ind_; }

    // Get the number of vectors that passed the test
    int howMany() const { return ind_.size(); }

    void setQuorum (int quorum) {
      state_ = Undefined;
      quorum_ = quorum;
    }

    int getQuorum() const { return quorum_; }

    void setTolerance(MagnitudeType tol)
    {
      state_ = Undefined;
      tol_ = tol;
    }

    MagnitudeType getTolerance() const { return tol_; }

    void setWhichNorm(ResType whichNorm)
    {
      state_ = Undefined;
      whichNorm_ = whichNorm;
    }

    ResType getWhichNorm() const { return whichNorm_; }

    void setScale(bool relscale)
    {
      state_ = Undefined;
      scaled_ = relscale;
    }

    bool getScale() const { return scaled_; }

    // Informs the status test that it should reset its internal configuration to the uninitialized state
    void reset()
    {
      ind_.resize(0);
      state_ = Undefined;
    }

    // Clears the results of the last status test
    void clearStatus() { reset(); };

    // Output formatted description of stopping test to output stream
    std::ostream & print(std::ostream &os, int indent=0) const;

  private:
    TestStatus state_;
    MagnitudeType tol_;
    std::vector<int> ind_;
    int quorum_;
    bool scaled_;
    enum ResType whichNorm_;
    bool throwExceptionOnNaN_;
    RCP<const OP> M_;

    const MagnitudeType ONE;
  };



  template <class ScalarType, class MV, class OP>
  StatusTestSpecTrans<ScalarType,MV,OP>::StatusTestSpecTrans(MagnitudeType tol, int quorum, ResType whichNorm, bool scaled, bool throwExceptionOnNaN, const RCP<const OP> Mop)
  : state_(Undefined),
    tol_(tol),
    quorum_(quorum),
    scaled_(scaled),
    whichNorm_(whichNorm),
    throwExceptionOnNaN_(throwExceptionOnNaN),
    M_(Mop),
    ONE(Teuchos::ScalarTraits<MagnitudeType>::one())
  {}



  template <class ScalarType, class MV, class OP>
  TestStatus StatusTestSpecTrans<ScalarType,MV,OP>::checkStatus( Eigensolver<ScalarType,MV,OP>* solver )
  {
    typedef Teuchos::ScalarTraits<MagnitudeType> MT;
    typedef TraceMinBase<ScalarType,MV,OP>       TS;

    // Cast the eigensolver to a TraceMin solver
    // Throw an exception if this fails
    TS* tm_solver = dynamic_cast<TS*>(solver);
    TEUCHOS_TEST_FOR_EXCEPTION(tm_solver == 0, std::invalid_argument, "The status test for spectral transformations currently only works for the trace minimization eigensolvers.  Sorry!");

    // Get the residual Ax-\lambda Bx, which is not computed when there's a spectral transformation
    // TraceMin computes Bx-1/\lambda Ax
    TraceMinBaseState<ScalarType,MV> state = tm_solver->getState();

    size_t nvecs = state.ritzShifts->size();
    std::vector<int> curind(nvecs);
    for(size_t i=0; i<nvecs; i++)
      curind[i] = i;

    RCP<const MV> locKX, locMX, locX;
    RCP<MV> R;
    locX = MVT::CloneView(*state.X,curind);
    if(state.KX != Teuchos::null)
      locKX = MVT::CloneView(*state.KX,curind);
    else
      locKX = locX;
    if(state.MX != Teuchos::null)
      locMX = MVT::CloneView(*state.MX,curind);
    else
      locMX = locX;
    R = MVT::CloneCopy(*locKX,curind);

    std::vector<MagnitudeType> evals(nvecs);
    for(size_t i=0; i<nvecs; i++)
      evals[i] = ONE/(*state.T)[i];
    MVT::MvScale(*R,evals);
    MVT::MvAddMv(-ONE,*R,ONE,*locMX,*R);

    // Compute the norms
    std::vector<MagnitudeType> res(nvecs);
    switch (whichNorm_) {
      case RES_2NORM:
      {
        MVT::MvNorm(*R,res);
        break;
      }
      case RES_ORTH:
      {
        RCP<MV> MR = MVT::Clone(*R,nvecs);
        OPT::Apply(*M_,*R,*MR);
        MVT::MvDot(*R,*MR,res);
        for(size_t i=0; i<nvecs; i++)
          res[i] = MT::squareroot(res[i]);
        break;
      }
      case RITZRES_2NORM:
      {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "The trace minimization solvers do not define a Ritz residual.  Please choose a different residual type");
        break;
      }
    }

    // if appropriate, scale the norms by the magnitude of the eigenvalue estimate
    if(scaled_)
    {
      for(size_t i=0; i<nvecs; i++)
        res[i] /= std::abs(evals[i]);
    }

    // test the norms
    ind_.resize(0);
    for(size_t i=0; i<nvecs; i++)
    {
      TEUCHOS_TEST_FOR_EXCEPTION( MT::isnaninf(res[i]), ResNormNaNError,
       "StatusTestSpecTrans::checkStatus(): residual norm is nan or inf" );
      if(res[i] < tol_)
        ind_.push_back(i);
    }
    int have = ind_.size();
    int need = (quorum_ == -1) ? nvecs : quorum_;
    state_ = (have >= need) ? Passed : Failed;
    return state_;
  }



  template <class ScalarType, class MV, class OP>
  std::ostream& StatusTestSpecTrans<ScalarType,MV,OP>::print(std::ostream& os, int indent) const
  {
    std::string ind(indent,' ');
    os << ind << "- StatusTestSpecTrans: ";
    switch (state_) {
      case Passed:
        os << "Passed\n";
        break;
      case Failed:
        os << "Failed\n";
        break;
      case Undefined:
        os << "Undefined\n";
        break;
    }
    os << ind << "  (Tolerance, WhichNorm,Scaled,Quorum): "
       << "(" << tol_;
    switch (whichNorm_) {
      case RES_ORTH:
        os << ",RES_ORTH";
        break;
      case RES_2NORM:
        os << ",RES_2NORM";
        break;
      case RITZRES_2NORM:
        os << ",RITZRES_2NORM";
        break;
    }
    os << "," << (scaled_ ? "true" : "false")
       << "," << quorum_
       << ")\n";

    if (state_ != Undefined) {
      os << ind << "  Which vectors: ";
      if (ind_.size() > 0) {
        for(size_t i=0; i<ind_.size(); i++) os << ind_[i] << " ";
        os << std::endl;
      }
      else
        os << "[empty]\n";
    }
    return os;
  }

}} // end of namespace

#endif
