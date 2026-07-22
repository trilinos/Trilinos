// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_VECTORPROFILER_H
#define ROL_VECTORPROFILER_H

#include "ROL_Vector.hpp"
#include <ostream>

namespace ROL {

/** @ingroup vector
    \class ROL::ProfiledVector 
    \brief By keeping a pointer to this in a derived 
           Vector class, a tally of all methods is kept
           for profiling function calls. 

    NOTE: This class is not yet compatible with vectors that have true duals

    In the cpp file where this is used, you must initialize the VectorFunctionCalls object.

    Example usage:

    template<>
    VectorFunctionCalls<int> ProfiledVector<int,double>::functionCalls = VectorFunctionCalls<int>();

*/

template<class Ordinal> 
struct VectorFunctionCalls {
  Ordinal constructor_;
  Ordinal destructor_;
  Ordinal plus_;
  Ordinal scale_;
  Ordinal dot_;
  Ordinal norm_;
  Ordinal clone_;
  Ordinal axpy_;
  Ordinal zero_;
  Ordinal basis_;
  Ordinal dimension_;
  Ordinal set_;
  Ordinal dual_; 
  Ordinal apply_;
  Ordinal applyUnary_;
  Ordinal applyBinary_;
  Ordinal reduce_;
  Ordinal setScalar_;
  Ordinal randomize_;
  VectorFunctionCalls() :
    constructor_(0), destructor_(0), plus_(0), scale_(0), dot_(0), norm_(0), clone_(0),
    axpy_(0), zero_(0), basis_(0), dimension_(0), set_(0), dual_(0), apply_(0),
    applyUnary_(0), applyBinary_(0), reduce_(0), setScalar_(0), randomize_(0) {}

}; // struct VectorFunctionCalls


// Forward declaration for friend functions
template<class Ordinal,class Real> 
class ProfiledVector;

template<class Ordinal,class Real>
VectorFunctionCalls<Ordinal> getVectorFunctionCalls( const ProfiledVector<Ordinal,Real> &x ) {
  return x.functionCalls_;
}

template<class Ordinal, class Real>
void printVectorFunctionCalls( const ProfiledVector<Ordinal,Real> &x, std::ostream &outStream = std::cout ) {
  outStream << "Total Vector Function Calls"  << std::endl;
  outStream << "---------------------------"  << std::endl;
  outStream << "Constructor : " << x.functionCalls_.constructor_  << std::endl;
  outStream << "Destructor  : " << x.functionCalls_.destructor_   << std::endl;
  outStream << "set         : " << x.functionCalls_.set_          << std::endl;
  outStream << "plus        : " << x.functionCalls_.plus_         << std::endl;
  outStream << "axpy        : " << x.functionCalls_.axpy_         << std::endl;
  outStream << "scale       : " << x.functionCalls_.scale_        << std::endl;
  outStream << "dot         : " << x.functionCalls_.dot_          << std::endl;
  outStream << "zero        : " << x.functionCalls_.zero_         << std::endl;
  outStream << "norm        : " << x.functionCalls_.norm_         << std::endl;
  outStream << "clone       : " << x.functionCalls_.clone_        << std::endl;
  outStream << "basis       : " << x.functionCalls_.basis_        << std::endl;
  outStream << "dual        : " << x.functionCalls_.dual_         << std::endl;
  outStream << "apply       : " << x.functionCalls_.apply_        << std::endl;
  outStream << "dimension   : " << x.functionCalls_.dimension_    << std::endl;
  outStream << "applyUnary  : " << x.functionCalls_.applyUnary_   << std::endl;
  outStream << "applyBinary : " << x.functionCalls_.applyBinary_  << std::endl;
  outStream << "reduce      : " << x.functionCalls_.reduce_       << std::endl;
  outStream << "setScalar   : " << x.functionCalls_.setScalar_    << std::endl;
  outStream << "randomize   : " << x.functionCalls_.randomize_    << std::endl;
}



template<class Ordinal,class Real> 
class ProfiledVector : public Vector<Real> {

  typedef Vector<Real>   V;

private: 
  ROL::Ptr<Vector<Real> > v_; 
  static VectorFunctionCalls<Ordinal> functionCalls_;
public:

  ProfiledVector( const ROL::Ptr<Vector<Real> > &v ) { 
    // Make sure that given vector is not itself a ProfiledVector to avoid recursion
    ROL::Ptr<ProfiledVector> pv = ROL::nullPtr;
    pv = ROL::dynamicPtrCast<ProfiledVector>(v);
    ROL_TEST_FOR_EXCEPTION( pv != ROL::nullPtr, std::logic_error, "ProfiledVector class "
    "cannot encapsulate a ProfiledVector object!");

    v_ = v;
    
    functionCalls_.constructor_++;
  }

  virtual ~ProfiledVector() {
    functionCalls_.destructor_++;
  }

  void plus( const Vector<Real> &x ) {
    ROL::Ptr<const V> xp = dynamic_cast<const ProfiledVector&>(x).getVector();

    functionCalls_.plus_++;
    v_->plus(*xp);
  } 

  void scale( const Real alpha ) {
    functionCalls_.scale_++;
    v_->scale(alpha);
  }
  
  Real dot( const Vector<Real> &x ) const {
    ROL::Ptr<const V> xp = dynamic_cast<const ProfiledVector&>(x).getVector();
    functionCalls_.dot_++;
    return v_->dot(*xp);
  }

  Real norm() const {
    functionCalls_.norm_++;
    return v_->norm();
  }

  ROL::Ptr<Vector<Real> > clone() const {
    functionCalls_.clone_++;
    return ROL::makePtr<ProfiledVector>( v_->clone() );
  }

  void axpy( const Real alpha, const Vector<Real> &x ) {
    ROL::Ptr<const V> xp = dynamic_cast<const ProfiledVector&>(x).getVector();
    functionCalls_.axpy_++;
    return v_->axpy(alpha,*xp);
  }

  void zero() {
    functionCalls_.zero_++;
    v_->zero();
  }

  ROL::Ptr<Vector<Real> > basis( const int i ) const {
    functionCalls_.basis_++;
   return ROL::makePtr<ProfiledVector>( v_->basis(i) );
  }

  int dimension() const {
    functionCalls_.dimension_++;
    return v_->dimension();
  }

  void set( const Vector<Real> &x ) {
    ROL::Ptr<const V> xp = dynamic_cast<const ProfiledVector&>(x).getVector();
    functionCalls_.set_++;
    v_->set(*xp);
  }

  // TODO: determine the correct way to handle dual when v_ is a generic ROL::Ptr<ROL::Vector>
  const Vector<Real> & dual() const {
    functionCalls_.dual_++; 
    return *this;
  }

  Real apply(const Vector<Real> &x) const {
    functionCalls_.apply_++;
    return v_->apply(x);
  }

  ROL::Ptr<Vector<Real> > getVector() {
    return v_;
  }

  ROL::Ptr<const Vector<Real> > getVector() const {
    return v_;
  }

  void applyUnary( const Elementwise::UnaryFunction<Real> &f ) {
    functionCalls_.applyUnary_++;
    v_->applyUnary(f);
  } 

  void applyBinary( const Elementwise::BinaryFunction<Real> &f, const Vector<Real> &x ) {
    functionCalls_.applyBinary_++;
    v_->applyBinary(f,x);
  }

  Real reduce( const Elementwise::ReductionOp<Real> &r ) const {
    functionCalls_.reduce_++;
    return v_->reduce(r);
  }

  void setScalar( const Real C ) {
    functionCalls_.setScalar_++;
    v_->setScalar(C);
  }

  void randomize( const Real l=0.0, const Real u=1.0) {
    functionCalls_.randomize_++;
    v_->randomize(l,u);
  }

  void print( std::ostream &outStream ) const {
    v_->print(outStream); 
  }
    
  friend VectorFunctionCalls<Ordinal> getVectorFunctionCalls<>( const ProfiledVector<Ordinal,Real> & );
  friend void printVectorFunctionCalls<>( const ProfiledVector<Ordinal,Real> &, std::ostream & );  

};


} // namespace ROL

#endif // ROL_RANDOMVECTOR_H
