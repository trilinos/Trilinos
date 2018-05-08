// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef ROL_RISKVECTOR_HPP
#define ROL_RISKVECTOR_HPP

#include "ROL_StdVector.hpp"
#include "ROL_RandVarFunctionalInfo.hpp"
#include "ROL_ParameterList.hpp"

namespace ROL {

template<class Real> 
class RiskVector : public Vector<Real> {
private:
  ROL::Ptr<std::vector<Real> > statObj_;
  ROL::Ptr<StdVector<Real> >   statObj_vec_;
  bool augmentedObj_;
  int nStatObj_;

  std::vector<ROL::Ptr<std::vector<Real> > > statCon_;
  std::vector<ROL::Ptr<StdVector<Real> > >   statCon_vec_;
  bool augmentedCon_;
  std::vector<int> nStatCon_;

  ROL::Ptr<Vector<Real> >      vec_;

  mutable bool isDualInitialized_;
  mutable ROL::Ptr<std::vector<Real> > dualObj_;
  mutable std::vector<ROL::Ptr<std::vector<Real> > > dualCon_;
  mutable ROL::Ptr<Vector<Real> > dual_vec1_;
  mutable ROL::Ptr<RiskVector<Real> > dual_vec_;

  void initializeObj(ROL::Ptr<ROL::ParameterList> &parlist,
               const Real stat = 1) {
    // Get risk measure information
    if (parlist != ROL::nullPtr) {
      std::string name;
      std::vector<Real> lower, upper;
      bool activated(false);
      RandVarFunctionalInfo<Real>(*parlist,name,nStatObj_,lower,upper,activated);
      augmentedObj_ = (nStatObj_ > 0) ? true : false;
      // Initialize statistic vector
      if (augmentedObj_) {
        statObj_     = ROL::makePtr<std::vector<Real>>(nStatObj_,stat);
        statObj_vec_ = ROL::makePtr<StdVector<Real>>(statObj_);
      }
    }
    else {
      augmentedObj_ = false;
      nStatObj_     = 0;
    }
  } 

  void initializeCon(std::vector<ROL::Ptr<ROL::ParameterList> > &parlist,
               const Real stat = 1) {
    int size = parlist.size();
    statCon_.resize(size); statCon_vec_.resize(size); nStatCon_.resize(size);
    for (int i = 0; i < size; ++i) { 
      if (parlist[i] != ROL::nullPtr) {
        // Get risk measure information
        std::string name;
        std::vector<Real> lower, upper;
        bool activated(false);
        RandVarFunctionalInfo<Real>(*parlist[i],name,nStatCon_[i],lower,upper,activated);
        augmentedCon_ = (nStatCon_[i] > 0) ? true : augmentedCon_;
        // Initialize statistic vector
        if (nStatCon_[i] > 0) {
          statCon_[i]     = ROL::makePtr<std::vector<Real>>(nStatCon_[i],stat);
          statCon_vec_[i] = ROL::makePtr<StdVector<Real>>(statCon_[i]);
        }
        else {
          statCon_[i]     = ROL::nullPtr;
          statCon_vec_[i] = ROL::nullPtr;
        }
      }
      else {
        statCon_[i]     = ROL::nullPtr;
        statCon_vec_[i] = ROL::nullPtr;
      }
    }
  }

public:
  
  // Objective risk only
  RiskVector( ROL::Ptr<ROL::ParameterList> &parlist,
        const ROL::Ptr<Vector<Real> >          &vec,
        const Real stat = 0 )
    : statObj_(ROL::nullPtr), statObj_vec_(ROL::nullPtr),
      augmentedObj_(false), nStatObj_(0),
      augmentedCon_(false),
      vec_(vec), isDualInitialized_(false) {
    initializeObj(parlist,stat);
  }

  // Inequality constraint risk only
  RiskVector( std::vector<ROL::Ptr<ROL::ParameterList> > &parlist,
        const ROL::Ptr<Vector<Real> > &vec,
        const Real stat = 0 )
    : statObj_(ROL::nullPtr), statObj_vec_(ROL::nullPtr),
      augmentedObj_(false), nStatObj_(0),
      augmentedCon_(false),
      vec_(vec), isDualInitialized_(false) {
    initializeCon(parlist,stat);
  }

  // Objective and inequality constraint risk
  RiskVector( ROL::Ptr<ROL::ParameterList> & parlistObj,
              std::vector<ROL::Ptr<ROL::ParameterList> > &parlistCon,
        const ROL::Ptr<Vector<Real> > &vec,
        const Real stat = 0 )
    : statObj_(ROL::nullPtr), statObj_vec_(ROL::nullPtr),
      augmentedObj_(false), nStatObj_(0),
      augmentedCon_(false),
      vec_(vec), isDualInitialized_(false) {
    initializeObj(parlistObj,stat);
    initializeCon(parlistCon,stat);
  }
 
  // Build from components
  RiskVector( const ROL::Ptr<Vector<Real> >                    &vec,
              const ROL::Ptr<std::vector<Real> >               &statObj,
              const std::vector<ROL::Ptr<std::vector<Real> > > &statCon )
    : statObj_(ROL::nullPtr), statObj_vec_(ROL::nullPtr),
      augmentedObj_(false), nStatObj_(0), augmentedCon_(false),
      vec_(vec), isDualInitialized_(false) {
    if (statObj != ROL::nullPtr) {
      statObj_      = statObj;
      statObj_vec_  = ROL::makePtr<StdVector<Real>>(statObj_);
      augmentedObj_ = true;
      nStatObj_     = statObj->size();
    }
    int size = statCon.size();
    statCon_.clear(); statCon_vec_.clear(); nStatCon_.clear();
    statCon_.resize(size,ROL::nullPtr);
    statCon_vec_.resize(size,ROL::nullPtr);
    nStatCon_.resize(size,0);
    for (int i = 0; i < size; ++i) {
      if (statCon[i] != ROL::nullPtr) {
        statCon_[i]     = statCon[i];
        statCon_vec_[i] = ROL::makePtr<StdVector<Real>>(statCon_[i]);
        augmentedCon_   = true;
        nStatCon_[i]    = statCon[i]->size();
      }
    }
  }

  // Build from components -- Objective only...no statistic
  RiskVector( const ROL::Ptr<Vector<Real> > &vec )
    : statObj_(ROL::nullPtr), statObj_vec_(ROL::nullPtr),
      augmentedObj_(false), nStatObj_(0), augmentedCon_(false),
      vec_(vec), isDualInitialized_(false) {}

  void set( const Vector<Real> &x ) {
    const RiskVector<Real> &xs = dynamic_cast<const RiskVector<Real>&>(x);
    vec_->set(*(xs.getVector()));
    if (augmentedObj_ && statObj_vec_ != ROL::nullPtr) {
      statObj_vec_->set(*(xs.getStatisticVector(0)));
    }
    if (augmentedCon_) {
      int size = statCon_vec_.size();
      for (int i = 0; i < size; ++i) {
        if (statCon_vec_[i] != ROL::nullPtr) {
          statCon_vec_[i]->set(*(xs.getStatisticVector(1,i)));
        }
      }
    }
  }

  void plus( const Vector<Real> &x ) {
    const RiskVector<Real> &xs = dynamic_cast<const RiskVector<Real>&>(x);
    vec_->plus(*(xs.getVector()));
    if (augmentedObj_ && statObj_vec_ != ROL::nullPtr) {
      statObj_vec_->plus(*(xs.getStatisticVector(0)));
    }
    if (augmentedCon_) {
      int size = statCon_vec_.size();
      for (int i = 0; i < size; ++i) {
        if (statCon_vec_[i] != ROL::nullPtr) {
          statCon_vec_[i]->plus(*(xs.getStatisticVector(1,i)));
        }
      }
    }
  }

  void scale( const Real alpha ) {
    vec_->scale(alpha);
    if (augmentedObj_ && statObj_vec_ != ROL::nullPtr) {
      statObj_vec_->scale(alpha);
    }
    if (augmentedCon_) {
      int size = statCon_vec_.size();
      for (int i = 0; i < size; ++i) {
        if (statCon_vec_[i] != ROL::nullPtr) {
          statCon_vec_[i]->scale(alpha);
        }
      }
    }
  }

  void axpy( const Real alpha, const Vector<Real> &x ) {
    const RiskVector<Real> &xs = dynamic_cast<const RiskVector<Real>&>(x);
    vec_->axpy(alpha,*(xs.getVector()));
    if (augmentedObj_ && statObj_vec_ != ROL::nullPtr) {
      statObj_vec_->axpy(alpha,*(xs.getStatisticVector(0)));
    }
    if (augmentedCon_) {
      int size = statCon_vec_.size();
      for (int i = 0; i < size; ++i) {
        if (statCon_vec_[i] != ROL::nullPtr) {
          statCon_vec_[i]->axpy(alpha,*(xs.getStatisticVector(1,i)));
        }
      }
    }
  }

  Real dot( const Vector<Real> &x ) const {
    const RiskVector<Real> &xs = dynamic_cast<const RiskVector<Real>&>(x);
    Real val = vec_->dot(*(xs.getVector()));
    if (augmentedObj_ && statObj_vec_ != ROL::nullPtr) {
      val += statObj_vec_->dot(*(xs.getStatisticVector(0)));
    }
    if (augmentedCon_) {
      int size = statCon_vec_.size();
      for (int i = 0; i < size; ++i) {
        if (statCon_vec_[i] != ROL::nullPtr) {
          val += statCon_vec_[i]->dot(*(xs.getStatisticVector(1,i)));
        }
      }
    }
    return val;
  }

  Real norm(void) const {
    return sqrt( dot(*this) );
  }

  ROL::Ptr<Vector<Real> > clone(void) const {
    ROL::Ptr<std::vector<Real> > e2 = ROL::nullPtr;
    if (augmentedObj_ && statObj_vec_ != ROL::nullPtr) {
      e2 = ROL::makePtr<std::vector<Real>>(nStatObj_,static_cast<Real>(0));
    }
    int size = statCon_vec_.size();
    std::vector<ROL::Ptr<std::vector<Real> > > e3(size, ROL::nullPtr);
    for (int j = 0; j < size; ++j) {
      if (statCon_vec_[j] != ROL::nullPtr) {
        e3[j] = ROL::makePtr<std::vector<Real>>(nStatCon_[j],static_cast<Real>(0));
      }
    }
    return ROL::makePtr<RiskVector>(vec_->clone(),e2,e3);
  }

  const Vector<Real> &dual(void) const {
    // Initialize dual vectors if not already initialized
    if ( !isDualInitialized_ ) {
      dual_vec1_ = vec_->dual().clone();
      dualObj_ = ROL::nullPtr;
      if (statObj_ != ROL::nullPtr) {
        dualObj_ = ROL::makePtr<std::vector<Real>>(statObj_->size());
      }
      int size = statCon_.size();
      dualCon_.clear(); dualCon_.resize(size,ROL::nullPtr);
      for (int i = 0; i < size; ++i) {
        if (statCon_[i] != ROL::nullPtr) {
          dualCon_[i] = ROL::makePtr<std::vector<Real>>(statCon_[i]->size());
        }
      }
      dual_vec_  = ROL::makePtr<RiskVector<Real>>(dual_vec1_,dualObj_,dualCon_);
      isDualInitialized_ = true;
    }
    // Set vector component 
    dual_vec1_->set(vec_->dual());
    // Set statistic component
    if ( augmentedObj_ && statObj_vec_ != ROL::nullPtr ) {
      ROL::dynamicPtrCast<RiskVector<Real> >(dual_vec_)->setStatistic(*statObj_,0);
    }
    if ( augmentedCon_ ) {
      int size = statCon_.size();
      for (int i = 0; i < size; ++i) {
        if (statCon_[i] != ROL::nullPtr) {
          ROL::dynamicPtrCast<RiskVector<Real> >(dual_vec_)->setStatistic(*statCon_[i],1,i);
        }
      }
    }
    // Return dual vector
    return *dual_vec_;
  }

  ROL::Ptr<Vector<Real> > basis( const int i )  const {
    ROL::Ptr<Vector<Real> > e1;
    ROL::Ptr<std::vector<Real> > e2 = ROL::nullPtr;
    if (augmentedObj_ && statObj_vec_ != ROL::nullPtr) {
      e2 = ROL::makePtr<std::vector<Real>>(nStatObj_,static_cast<Real>(0));
    }
    int size = statCon_vec_.size();
    std::vector<ROL::Ptr<std::vector<Real> > > e3(size);
    for (int j = 0; j < size; ++j) {
      if (statCon_vec_[j] != ROL::nullPtr) {
        e3[j] = ROL::makePtr<std::vector<Real>>(nStatCon_[j],static_cast<Real>(0));
      }
    }
    int n1 = vec_->dimension(), n2 = 0;
    if (statObj_vec_ != ROL::nullPtr) {
      n2 = statObj_vec_->dimension();
    }
    if ( i < n1 ) {
      e1 = vec_->basis(i);
    }
    else if (i >= n1 && i < n1+n2) {
      e1 = vec_->clone(); e1->zero();
      (*e2)[i-n1] = static_cast<Real>(1);
    }
    else if (i >= n1+n2) {
      e1 = vec_->clone(); e1->zero();
      int sum = n1+n2, sum0 = sum;
      for (int j = 0; j < size; ++j) {
        if (statCon_vec_[j] != ROL::nullPtr) {
          sum += nStatCon_[j];
          if (i < sum) {
            (*e3[j])[i-sum0] = static_cast<Real>(1);
            break;
          }
          sum0 = sum;
        }
      }
      if (i >= sum) {
        throw Exception::NotImplemented(">>> ROL::RiskVector::Basis: index out of bounds!");
      }
    }
    return ROL::makePtr<RiskVector<Real>>(e1,e2,e3);
  }

  void applyUnary( const Elementwise::UnaryFunction<Real> &f ) {
    vec_->applyUnary(f);
    if (augmentedObj_ && statObj_vec_ != ROL::nullPtr) {
      statObj_vec_->applyUnary(f);
    }
    if (augmentedCon_) {
      int size = statCon_vec_.size();
      for (int i = 0; i < size; ++i) {
        if (statCon_vec_[i] != ROL::nullPtr) {
          statCon_vec_[i]->applyUnary(f);
        }
      }
    }
  }

  void applyBinary( const Elementwise::BinaryFunction<Real> &f, const Vector<Real> &x ) {
    const RiskVector<Real> &xs = dynamic_cast<const RiskVector<Real>&>(x);
    vec_->applyBinary(f,*xs.getVector());
    if (augmentedObj_ && statObj_vec_ != ROL::nullPtr) {
      statObj_vec_->applyBinary(f,*xs.getStatisticVector(0));
    }
    if (augmentedCon_) {
      int size = statCon_vec_.size();
      for (int i = 0; i < size; ++i) {
        if (statCon_vec_[i] != ROL::nullPtr) {
          statCon_vec_[i]->applyBinary(f,*xs.getStatisticVector(1,i));
        }
      }
    }
  }

  Real reduce( const Elementwise::ReductionOp<Real> &r ) const {
    Real result = r.initialValue();
    r.reduce(vec_->reduce(r),result);
    if (augmentedObj_ && statObj_vec_ != ROL::nullPtr) {
      r.reduce(statObj_vec_->reduce(r),result);
    }
    if (augmentedCon_) {
      int size = statCon_vec_.size();
      for (int i = 0; i < size; ++i) {
        if (statCon_vec_[i] != ROL::nullPtr) {
          r.reduce(statCon_vec_[i]->reduce(r),result);
        }
      }
    }
    return result;
  }

  void setScalar( const Real C ) {
    vec_->setScalar(C);
    if (augmentedObj_ && statObj_vec_ != ROL::nullPtr) {
      statObj_vec_->setScalar(C);
    }
    if (augmentedCon_) {
      int size = statCon_vec_.size();
      for (int i = 0; i < size; ++i) {
        if (statCon_vec_[i] != ROL::nullPtr) {
          statCon_vec_[i]->setScalar(C);
        }
      }
    }
  }

  int dimension(void) const {
    int dim = vec_->dimension();
    if (augmentedObj_) {
      dim += statObj_vec_->dimension();
    }
    if (augmentedCon_) {
      int size = statCon_vec_.size();
      for (int i = 0; i < size; ++i) {
        if (statCon_vec_[i] != ROL::nullPtr) {
          dim += statCon_vec_[i]->dimension();
        }
      }
    }
    return dim;
  }

  /***************************************************************************/
  /************ ROL VECTOR ACCESSOR FUNCTIONS ********************************/
  /***************************************************************************/
  ROL::Ptr<const StdVector<Real>> 
  getStatisticVector(const int comp, const int index = 0) const {
    if (comp == 0) {
      return statObj_vec_;
    }
    else if (comp == 1) {
      return statCon_vec_[index];
    }
    else {
      throw Exception::NotImplemented(">>> ROL::RiskVector::getStatisticVector: Component must be 0 or 1!");
    }
    return ROL::nullPtr;
  }

  ROL::Ptr<StdVector<Real>> 
  getStatisticVector(const int comp, const int index = 0) {
    if (comp == 0) {
      return statObj_vec_;
    }
    else if (comp == 1) {
      return statCon_vec_[index];
    }
    else {
      throw Exception::NotImplemented(">>> ROL::RiskVector::getStatistic: Component must be 0 or 1!");
    }
    return ROL::nullPtr;
  }

  ROL::Ptr<const Vector<Real> > getVector(void) const {
    return vec_;
  }

  ROL::Ptr<Vector<Real> > getVector(void) {
    return vec_;
  }

  /***************************************************************************/
  /************ COMPONENT ACCESSOR FUNCTIONS *********************************/
  /***************************************************************************/
  ROL::Ptr<std::vector<Real>> 
  getStatistic(const int comp = 0, const int index = 0) {
    if (comp == 0) {
      if (augmentedObj_) {
        return statObj_;
      }
    }
    else if (comp == 1) {
      if (augmentedCon_) {
        return statCon_[index];
      }
    }
    else {
      throw Exception::NotImplemented(">>> ROL::RiskVector::getStatistic: Component must be 0 or 1!");
    }
    return ROL::nullPtr;
  }

  ROL::Ptr<const std::vector<Real>> 
  getStatistic(const int comp = 0, const int index = 0) const {
    if (comp == 0) {
      if (augmentedObj_) {
        return statObj_;
      }
    }
    else if (comp == 1) {
      if (augmentedCon_) {
        return statCon_[index];
      }
    }
    else {
      throw Exception::NotImplemented(">>> ROL::RiskVector::getStatistic: Component must be 0 or 1!");
    }
    return ROL::nullPtr;
  }

  void setStatistic(const Real stat, const int comp = 0, const int index = 0) {
    if ( comp == 0 ) {
      if ( augmentedObj_ ) {
        statObj_->assign(nStatObj_,stat);
      }
    }
    else if ( comp == 1 ) {
      if ( augmentedCon_ ) {
        statCon_[index]->assign(nStatCon_[index],stat);
      }
    }
    else {
      throw Exception::NotImplemented(">>> ROL::RiskVector::setStatistic: Component must be 0 or 1!");
    }
  }
 
  void setStatistic(const std::vector<Real> &stat, const int comp = 0, const int index = 0) {
    if ( comp == 0 ) {
      if ( augmentedObj_ ) {
        if ( nStatObj_ != static_cast<int>(stat.size()) ) {
          throw Exception::NotImplemented(">>> ROL::RiskVector::setStatistic: Dimension mismatch!");
        }
        statObj_->assign(stat.begin(),stat.end());
      }
    }
    else if ( comp == 1) {
      if ( augmentedCon_ ) {
        if ( nStatCon_[index] != static_cast<int>(stat.size()) ) {
          throw Exception::NotImplemented(">>> ROL::RiskVector::setStatistic: Dimension mismatch!");
        }
        statCon_[index]->assign(stat.begin(),stat.end());
      }
    }
    else {
      throw Exception::NotImplemented(">>> ROL::RiskVector::setStatistic: Component must be 0 or 1!");
    }
  }
 
  void setVector(const Vector<Real>& vec) {
    vec_->set(vec);
  }
};

}

#endif
