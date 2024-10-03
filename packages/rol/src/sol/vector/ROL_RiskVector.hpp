// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
  Ptr<std::vector<Real> > statObj_;
  Ptr<StdVector<Real> >   statObj_vec_;
  bool augmentedObj_;
  int nStatObj_;

  std::vector<Ptr<std::vector<Real> > > statCon_;
  std::vector<Ptr<StdVector<Real> > >   statCon_vec_;
  bool augmentedCon_;
  std::vector<int> nStatCon_;

  Ptr<Vector<Real> >      vec_;

  mutable bool isDualInitialized_;
  mutable Ptr<std::vector<Real> > dualObj_;
  mutable std::vector<Ptr<std::vector<Real> > > dualCon_;
  mutable Ptr<Vector<Real> > dual_vec1_;
  mutable Ptr<RiskVector<Real> > dual_vec_;

  void initializeObj(Ptr<ParameterList> &parlist,
               const Real stat = 1) {
    // Get risk measure information
    if (parlist != nullPtr) {
      std::string name;
      std::vector<Real> lower, upper;
      bool activated(false);
      RandVarFunctionalInfo<Real>(*parlist,name,nStatObj_,lower,upper,activated);
      augmentedObj_ = (nStatObj_ > 0) ? true : false;
      // Initialize statistic vector
      if (augmentedObj_) {
        statObj_     = makePtr<std::vector<Real>>(nStatObj_,stat);
        statObj_vec_ = makePtr<StdVector<Real>>(statObj_);
      }
    }
    else {
      augmentedObj_ = false;
      nStatObj_     = 0;
    }
  }

  void initializeCon(std::vector<Ptr<ParameterList> > &parlist,
               const Real stat = 1) {
    int size = parlist.size();
    statCon_.resize(size); statCon_vec_.resize(size); nStatCon_.resize(size);
    for (int i = 0; i < size; ++i) { 
      if (parlist[i] != nullPtr) {
        // Get risk measure information
        std::string name;
        std::vector<Real> lower, upper;
        bool activated(false);
        RandVarFunctionalInfo<Real>(*parlist[i],name,nStatCon_[i],lower,upper,activated);
        augmentedCon_ = (nStatCon_[i] > 0) ? true : augmentedCon_;
        // Initialize statistic vector
        if (nStatCon_[i] > 0) {
          statCon_[i]     = makePtr<std::vector<Real>>(nStatCon_[i],stat);
          statCon_vec_[i] = makePtr<StdVector<Real>>(statCon_[i]);
        }
        else {
          statCon_[i]     = nullPtr;
          statCon_vec_[i] = nullPtr;
        }
      }
      else {
        statCon_[i]     = nullPtr;
        statCon_vec_[i] = nullPtr;
      }
    }
  }

public:
  
  // Objective risk only
  RiskVector( Ptr<ParameterList> &parlist,
        const Ptr<Vector<Real> >          &vec,
        const Real stat = 0 )
    : statObj_(nullPtr), statObj_vec_(nullPtr),
      augmentedObj_(false), nStatObj_(0),
      augmentedCon_(false),
      vec_(vec), isDualInitialized_(false) {
    initializeObj(parlist,stat);
  }

  // Inequality constraint risk only
  RiskVector( std::vector<Ptr<ParameterList> > &parlist,
        const Ptr<Vector<Real> > &vec,
        const Real stat = 0 )
    : statObj_(nullPtr), statObj_vec_(nullPtr),
      augmentedObj_(false), nStatObj_(0),
      augmentedCon_(false),
      vec_(vec), isDualInitialized_(false) {
    initializeCon(parlist,stat);
  }

  // Objective and inequality constraint risk
  RiskVector( Ptr<ParameterList> & parlistObj,
              std::vector<Ptr<ParameterList> > &parlistCon,
        const Ptr<Vector<Real> > &vec,
        const Real stat = 0 )
    : statObj_(nullPtr), statObj_vec_(nullPtr),
      augmentedObj_(false), nStatObj_(0),
      augmentedCon_(false),
      vec_(vec), isDualInitialized_(false) {
    initializeObj(parlistObj,stat);
    initializeCon(parlistCon,stat);
  }
 
  // Build from components
  RiskVector( const Ptr<Vector<Real> >                    &vec,
              const Ptr<std::vector<Real> >               &statObj,
              const std::vector<Ptr<std::vector<Real> > > &statCon )
    : statObj_(nullPtr), statObj_vec_(nullPtr),
      augmentedObj_(false), nStatObj_(0), augmentedCon_(false),
      vec_(vec), isDualInitialized_(false) {
    if (statObj != nullPtr) {
      statObj_      = statObj;
      statObj_vec_  = makePtr<StdVector<Real>>(statObj_);
      augmentedObj_ = true;
      nStatObj_     = statObj->size();
    }
    int size = statCon.size();
    statCon_.clear(); statCon_vec_.clear(); nStatCon_.clear();
    statCon_.resize(size,nullPtr);
    statCon_vec_.resize(size,nullPtr);
    nStatCon_.resize(size,0);
    for (int i = 0; i < size; ++i) {
      if (statCon[i] != nullPtr) {
        statCon_[i]     = statCon[i];
        statCon_vec_[i] = makePtr<StdVector<Real>>(statCon_[i]);
        augmentedCon_   = true;
        nStatCon_[i]    = statCon[i]->size();
      }
    }
  }

  // Build from components -- Objective only...no statistic
  RiskVector( const Ptr<Vector<Real> > &vec )
    : statObj_(nullPtr), statObj_vec_(nullPtr),
      augmentedObj_(false), nStatObj_(0), augmentedCon_(false),
      vec_(vec), isDualInitialized_(false) {}

  void set( const Vector<Real> &x ) {
    const RiskVector<Real> &xs = dynamic_cast<const RiskVector<Real>&>(x);
    vec_->set(*(xs.getVector()));
    if (augmentedObj_ && statObj_vec_ != nullPtr) {
      statObj_vec_->set(*(xs.getStatisticVector(0)));
    }
    if (augmentedCon_) {
      int size = statCon_vec_.size();
      for (int i = 0; i < size; ++i) {
        if (statCon_vec_[i] != nullPtr) {
          statCon_vec_[i]->set(*(xs.getStatisticVector(1,i)));
        }
      }
    }
  }

  void plus( const Vector<Real> &x ) {
    const RiskVector<Real> &xs = dynamic_cast<const RiskVector<Real>&>(x);
    vec_->plus(*(xs.getVector()));
    if (augmentedObj_ && statObj_vec_ != nullPtr) {
      statObj_vec_->plus(*(xs.getStatisticVector(0)));
    }
    if (augmentedCon_) {
      int size = statCon_vec_.size();
      for (int i = 0; i < size; ++i) {
        if (statCon_vec_[i] != nullPtr) {
          statCon_vec_[i]->plus(*(xs.getStatisticVector(1,i)));
        }
      }
    }
  }

  void scale( const Real alpha ) {
    vec_->scale(alpha);
    if (augmentedObj_ && statObj_vec_ != nullPtr) {
      statObj_vec_->scale(alpha);
    }
    if (augmentedCon_) {
      int size = statCon_vec_.size();
      for (int i = 0; i < size; ++i) {
        if (statCon_vec_[i] != nullPtr) {
          statCon_vec_[i]->scale(alpha);
        }
      }
    }
  }

  void axpy( const Real alpha, const Vector<Real> &x ) {
    const RiskVector<Real> &xs = dynamic_cast<const RiskVector<Real>&>(x);
    vec_->axpy(alpha,*(xs.getVector()));
    if (augmentedObj_ && statObj_vec_ != nullPtr) {
      statObj_vec_->axpy(alpha,*(xs.getStatisticVector(0)));
    }
    if (augmentedCon_) {
      int size = statCon_vec_.size();
      for (int i = 0; i < size; ++i) {
        if (statCon_vec_[i] != nullPtr) {
          statCon_vec_[i]->axpy(alpha,*(xs.getStatisticVector(1,i)));
        }
      }
    }
  }

  Real dot( const Vector<Real> &x ) const {
    const RiskVector<Real> &xs = dynamic_cast<const RiskVector<Real>&>(x);
    Real val = vec_->dot(*(xs.getVector()));
    if (augmentedObj_ && statObj_vec_ != nullPtr) {
      val += statObj_vec_->dot(*(xs.getStatisticVector(0)));
    }
    if (augmentedCon_) {
      int size = statCon_vec_.size();
      for (int i = 0; i < size; ++i) {
        if (statCon_vec_[i] != nullPtr) {
          val += statCon_vec_[i]->dot(*(xs.getStatisticVector(1,i)));
        }
      }
    }
    return val;
  }

  Real norm(void) const {
    return sqrt( dot(*this) );
  }

  Ptr<Vector<Real> > clone(void) const {
    Ptr<std::vector<Real> > e2 = nullPtr;
    if (augmentedObj_ && statObj_vec_ != nullPtr) {
      e2 = makePtr<std::vector<Real>>(nStatObj_,static_cast<Real>(0));
    }
    int size = statCon_vec_.size();
    std::vector<Ptr<std::vector<Real> > > e3(size, nullPtr);
    for (int j = 0; j < size; ++j) {
      if (statCon_vec_[j] != nullPtr) {
        e3[j] = makePtr<std::vector<Real>>(nStatCon_[j],static_cast<Real>(0));
      }
    }
    return makePtr<RiskVector>(vec_->clone(),e2,e3);
  }

  const Vector<Real> &dual(void) const {
    // Initialize dual vectors if not already initialized
    if ( !isDualInitialized_ ) {
      dual_vec1_ = vec_->dual().clone();
      dualObj_ = nullPtr;
      if (statObj_ != nullPtr) {
        dualObj_ = makePtr<std::vector<Real>>(statObj_->size());
      }
      int size = statCon_.size();
      dualCon_.clear(); dualCon_.resize(size,nullPtr);
      for (int i = 0; i < size; ++i) {
        if (statCon_[i] != nullPtr) {
          dualCon_[i] = makePtr<std::vector<Real>>(statCon_[i]->size());
        }
      }
      dual_vec_  = makePtr<RiskVector<Real>>(dual_vec1_,dualObj_,dualCon_);
      isDualInitialized_ = true;
    }
    // Set vector component 
    dual_vec1_->set(vec_->dual());
    // Set statistic component
    if ( augmentedObj_ && statObj_vec_ != nullPtr ) {
      dynamicPtrCast<RiskVector<Real> >(dual_vec_)->setStatistic(*statObj_,0);
    }
    if ( augmentedCon_ ) {
      int size = statCon_.size();
      for (int i = 0; i < size; ++i) {
        if (statCon_[i] != nullPtr) {
          dynamicPtrCast<RiskVector<Real> >(dual_vec_)->setStatistic(*statCon_[i],1,i);
        }
      }
    }
    // Return dual vector
    return *dual_vec_;
  }

  Ptr<Vector<Real> > basis( const int i )  const {
    Ptr<Vector<Real> > e1;
    Ptr<std::vector<Real> > e2 = nullPtr;
    if (augmentedObj_ && statObj_vec_ != nullPtr) {
      e2 = makePtr<std::vector<Real>>(nStatObj_,static_cast<Real>(0));
    }
    int size = statCon_vec_.size();
    std::vector<Ptr<std::vector<Real> > > e3(size);
    for (int j = 0; j < size; ++j) {
      if (statCon_vec_[j] != nullPtr) {
        e3[j] = makePtr<std::vector<Real>>(nStatCon_[j],static_cast<Real>(0));
      }
    }
    int n1 = vec_->dimension(), n2 = 0;
    if (statObj_vec_ != nullPtr) {
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
        if (statCon_vec_[j] != nullPtr) {
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
    return makePtr<RiskVector<Real>>(e1,e2,e3);
  }

  void applyUnary( const Elementwise::UnaryFunction<Real> &f ) {
    vec_->applyUnary(f);
    if (augmentedObj_ && statObj_vec_ != nullPtr) {
      statObj_vec_->applyUnary(f);
    }
    if (augmentedCon_) {
      int size = statCon_vec_.size();
      for (int i = 0; i < size; ++i) {
        if (statCon_vec_[i] != nullPtr) {
          statCon_vec_[i]->applyUnary(f);
        }
      }
    }
  }

  void applyBinary( const Elementwise::BinaryFunction<Real> &f, const Vector<Real> &x ) {
    const RiskVector<Real> &xs = dynamic_cast<const RiskVector<Real>&>(x);
    vec_->applyBinary(f,*xs.getVector());
    if (augmentedObj_ && statObj_vec_ != nullPtr) {
      statObj_vec_->applyBinary(f,*xs.getStatisticVector(0));
    }
    if (augmentedCon_) {
      int size = statCon_vec_.size();
      for (int i = 0; i < size; ++i) {
        if (statCon_vec_[i] != nullPtr) {
          statCon_vec_[i]->applyBinary(f,*xs.getStatisticVector(1,i));
        }
      }
    }
  }

  Real reduce( const Elementwise::ReductionOp<Real> &r ) const {
    Real result = r.initialValue();
    r.reduce(vec_->reduce(r),result);
    if (augmentedObj_ && statObj_vec_ != nullPtr) {
      r.reduce(statObj_vec_->reduce(r),result);
    }
    if (augmentedCon_) {
      int size = statCon_vec_.size();
      for (int i = 0; i < size; ++i) {
        if (statCon_vec_[i] != nullPtr) {
          r.reduce(statCon_vec_[i]->reduce(r),result);
        }
      }
    }
    return result;
  }

  void setScalar( const Real C ) {
    vec_->setScalar(C);
    if (augmentedObj_ && statObj_vec_ != nullPtr) {
      statObj_vec_->setScalar(C);
    }
    if (augmentedCon_) {
      int size = statCon_vec_.size();
      for (int i = 0; i < size; ++i) {
        if (statCon_vec_[i] != nullPtr) {
          statCon_vec_[i]->setScalar(C);
        }
      }
    }
  }

  void randomize( const Real l=0.0, const Real u=1.0 ) {
    vec_->randomize(l,u);
    if (augmentedObj_ && statObj_vec_ != nullPtr) {
      statObj_vec_->randomize(l,u);
    }
    if (augmentedCon_) {
      int size = statCon_vec_.size();
      for (int i = 0; i < size; ++i) {
        if (statCon_vec_[i] != nullPtr) {
          statCon_vec_[i]->randomize(l,u);
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
        if (statCon_vec_[i] != nullPtr) {
          dim += statCon_vec_[i]->dimension();
        }
      }
    }
    return dim;
  }

  /***************************************************************************/
  /************ ROL VECTOR ACCESSOR FUNCTIONS ********************************/
  /***************************************************************************/
  Ptr<const StdVector<Real>> 
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
  }

  Ptr<StdVector<Real>> 
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
  }

  Ptr<const Vector<Real> > getVector(void) const {
    return vec_;
  }

  Ptr<Vector<Real> > getVector(void) {
    return vec_;
  }

  /***************************************************************************/
  /************ COMPONENT ACCESSOR FUNCTIONS *********************************/
  /***************************************************************************/
  Ptr<std::vector<Real>> 
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
    return nullPtr;
  }

  Ptr<const std::vector<Real>> 
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
    return nullPtr;
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
