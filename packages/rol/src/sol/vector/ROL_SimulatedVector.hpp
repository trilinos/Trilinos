// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "ROL_Vector.hpp"
#include "ROL_SampleGenerator.hpp"

#ifndef ROL_SIMULATED_VECTOR_H
#define ROL_SIMULATED_VECTOR_H

/** @ingroup la_group
 *  \class ROL::SimulatedVector
 *  \brief Defines the linear algebra of a vector space on a generic
           partitioned vector where the individual vectors are
           distributed in batches defined by ROL::BatchManager.
           This is a batch-distributed version of ROL::PartitionedVector.
 */


namespace ROL {

template<class Real>
class PrimalSimulatedVector;

template<class Real>
class DualSimulatedVector;

template<class Real>
class SimulatedVector : public Vector<Real> {

  typedef Vector<Real>                   V;
  typedef ROL::Ptr<V>                    Vp;
  typedef ROL::Ptr<BatchManager<Real> >  VBMp; 
  typedef SimulatedVector<Real>          PV;

private:
  const std::vector<Vp>                  vecs_;
  ROL::Ptr<BatchManager<Real> >          bman_;
  mutable std::vector<Vp>                dual_vecs_;
  mutable ROL::Ptr<PV>                   dual_pvec_;
public:

  typedef typename std::vector<PV>::size_type    size_type;

  SimulatedVector( const std::vector<Vp> &vecs, const VBMp &bman ) : vecs_(vecs), bman_(bman) {
    for( size_type i=0; i<vecs_.size(); ++i ) {
      dual_vecs_.push_back((vecs_[i]->dual()).clone());
    }
  }

  void set( const V &x ) {
    
    const PV &xs = dynamic_cast<const PV&>(x);

    ROL_TEST_FOR_EXCEPTION( numVectors() != xs.numVectors(),
                            std::invalid_argument,
                            "Error: Vectors must have the same number of subvectors." );

    for( size_type i=0; i<vecs_.size(); ++i ) {
      vecs_[i]->set(*xs.get(i));
    }
  }

  void plus( const V &x ) {
    
    const PV &xs = dynamic_cast<const PV&>(x);

    ROL_TEST_FOR_EXCEPTION( numVectors() != xs.numVectors(),
                            std::invalid_argument,
                            "Error: Vectors must have the same number of subvectors." );

    for( size_type i=0; i<vecs_.size(); ++i ) {
      vecs_[i]->plus(*xs.get(i));
    }
  }

  void scale( const Real alpha ) {
    for( size_type i=0; i<vecs_.size(); ++i ) {
      vecs_[i]->scale(alpha);
    }
  }

  void axpy( const Real alpha, const V &x ) {
    
    const PV &xs = dynamic_cast<const PV&>(x);

    ROL_TEST_FOR_EXCEPTION( numVectors() != xs.numVectors(),
                            std::invalid_argument,
                            "Error: Vectors must have the same number of subvectors." );

    for( size_type i=0; i<vecs_.size(); ++i ) {
      vecs_[i]->axpy(alpha,*xs.get(i));
    }
  }

  virtual Real dot( const V &x ) const {
    
    const PV &xs = dynamic_cast<const PV&>(x);

   ROL_TEST_FOR_EXCEPTION( numVectors() != xs.numVectors(),
                           std::invalid_argument,
                           "Error: Vectors must have the same number of subvectors." );

    Real locresult = 0;
    Real result = 0;
    for( size_type i=0; i<vecs_.size(); ++i ) {
      locresult += vecs_[i]->dot(*xs.get(i));
    }

    bman_->sumAll(&locresult, &result, 1);

    return result;
  }

  Real norm() const {
    return std::sqrt(dot(*this));
  }

  virtual Vp clone() const {
    
    

    std::vector<Vp> clonevec;
    for( size_type i=0; i<vecs_.size(); ++i ) {
      clonevec.push_back(vecs_[i]->clone());
    }
    return ROL::makePtr<PV>(clonevec, bman_);
  }

  virtual const V& dual(void) const {
    

    for( size_type i=0; i<vecs_.size(); ++i ) {
      dual_vecs_[i]->set(vecs_[i]->dual());
    }
    dual_pvec_ = ROL::makePtr<PV>( dual_vecs_, bman_ );
    return *dual_pvec_;
  }

  Vp basis( const int i ) const { // this must be fixed for distributed batching

    ROL_TEST_FOR_EXCEPTION( i >= dimension() || i<0,
                            std::invalid_argument,
                            "Error: Basis index must be between 0 and vector dimension." );
    Vp bvec = clone();

    // Downcast
    PV &eb = dynamic_cast<PV&>(*bvec);

    int begin = 0;
    int end = 0;

    // Iterate over subvectors
    for( size_type j=0; j<vecs_.size(); ++j ) {

      end += vecs_[j]->dimension();

      if( begin<= i && i<end ) {
        eb.set(j, *(vecs_[j]->basis(i-begin)) );
      }
      else {
        eb.zero(j);
      }

      begin = end;

    }
    return bvec;
  }

  int dimension() const { // this must be fixed for distributed batching
    int total_dim = 0;
    for( size_type j=0; j<vecs_.size(); ++j ) {
      total_dim += vecs_[j]->dimension();
    }
    return total_dim;
  }

  void zero() {
    for( size_type j=0; j<vecs_.size(); ++j ) {
      vecs_[j]->zero();
    }
  }

  // Apply the same unary function to each subvector
  void applyUnary( const Elementwise::UnaryFunction<Real> &f ) {
    for( size_type i=0; i<vecs_.size(); ++i ) {
      vecs_[i]->applyUnary(f);
    }
  }

  // Apply the same binary function to each pair of subvectors in this vector and x
  void applyBinary( const Elementwise::BinaryFunction<Real> &f, const V &x ) {
    const PV &xs = dynamic_cast<const PV&>(x);

    for( size_type i=0; i<vecs_.size(); ++i ) {
      vecs_[i]->applyBinary(f,*xs.get(i));
    }
  }

  Real reduce( const Elementwise::ReductionOp<Real> &r ) const {
    Real result = r.initialValue();

    for( size_type i=0; i<vecs_.size(); ++i ) {
      r.reduce(vecs_[i]->reduce(r),result);
    }
    return result;
  }

  void setScalar(const Real C) {
    for( size_type i=0; i<vecs_.size(); ++i ) {
      vecs_[i]->setScalar(C);
    }
  }

  void randomize(const Real l=0.0, const Real u=1.0) {
    for( size_type i=0; i<vecs_.size(); ++i ) {
      vecs_[i]->randomize(l,u);
    }
  }

  // Methods that do not exist in the base class

  // In distributed batching mode, these are understood to take local indices.

  ROL::Ptr<const Vector<Real> > get(size_type i) const {
    return vecs_[i];
  }

  ROL::Ptr<Vector<Real> > get(size_type i) {
    return vecs_[i];
  }

  void set(size_type i, const V &x) {
    vecs_[i]->set(x);
  }

  void zero(size_type i) {
    vecs_[i]->zero();
  }

  size_type numVectors() const {
    return vecs_.size();
  }

};

// Helper methods
template<class Real>
ROL::Ptr<Vector<Real> > 
CreateSimulatedVector( const ROL::Ptr<Vector<Real>> &a, 
                       const ROL::Ptr<BatchManager<Real>> &bman ) {
  
  typedef ROL::Ptr<Vector<Real> >       Vp;
  typedef SimulatedVector<Real>  PV;

  Vp temp[] = {a};
  return ROL::makePtr<PV>(std::vector<Vp>(temp, temp+1), bman );
}

template<class Real>
class PrimalSimulatedVector : public SimulatedVector<Real> {
private:
  const std::vector<ROL::Ptr<Vector<Real>>>    vecs_;
  const ROL::Ptr<BatchManager<Real>>           bman_;
  const ROL::Ptr<SampleGenerator<Real>>        sampler_;
  mutable std::vector<ROL::Ptr<Vector<Real>>>  dual_vecs_;
  mutable ROL::Ptr<DualSimulatedVector<Real>>  dual_pvec_;
  mutable bool isDualInitialized_;
public:

  PrimalSimulatedVector(const std::vector<ROL::Ptr<Vector<Real>>> &vecs,
                        const ROL::Ptr<BatchManager<Real>>         &bman,
                        const ROL::Ptr<SampleGenerator<Real>>      &sampler)
    : SimulatedVector<Real>(vecs,bman), vecs_(vecs), bman_(bman), sampler_(sampler),
      isDualInitialized_(false) {
    for( int i=0; i<sampler_->numMySamples(); ++i ) {
      dual_vecs_.push_back((vecs_[i]->dual()).clone());
    }
  }

  Real dot(const Vector<Real> &x) const {
    const SimulatedVector<Real> &xs
      = dynamic_cast<const SimulatedVector<Real>&>(x);

   ROL_TEST_FOR_EXCEPTION( sampler_->numMySamples() != static_cast<int>(xs.numVectors()),
                               std::invalid_argument,
                               "Error: Vectors must have the same number of subvectors." );

    Real c = 0;
    Real locresult = 0;
    Real result = 0;
    for( int i=0; i<sampler_->numMySamples(); ++i ) {
      //locresult += sampler_->getMyWeight(i) * vecs_[i]->dot(*xs.get(i));
      Real y = sampler_->getMyWeight(i) * vecs_[i]->dot(*xs.get(i)) - c;
      Real t = locresult + y;
      c = (t - locresult) - y;
      locresult = t;
    }

    bman_->sumAll(&locresult, &result, 1);

    return result;
  }

  ROL::Ptr<Vector<Real> > clone(void) const {
    std::vector<ROL::Ptr<Vector<Real> > > clonevec;
    for( int i=0; i<sampler_->numMySamples(); ++i ) {
      clonevec.push_back(vecs_[i]->clone());
    }
    return ROL::makePtr<PrimalSimulatedVector<Real>>(clonevec, bman_, sampler_);
  }

  const Vector<Real>& dual(void) const {
    if (!isDualInitialized_) {
      dual_pvec_ = ROL::makePtr<DualSimulatedVector<Real>>(dual_vecs_, bman_, sampler_);
      isDualInitialized_ = true;
    }
    for( int i=0; i<sampler_->numMySamples(); ++i ) {
      dual_vecs_[i]->set(vecs_[i]->dual());
      dual_vecs_[i]->scale(sampler_->getMyWeight(i));
    }
    return *dual_pvec_;
  }

};

template<class Real>
class DualSimulatedVector : public SimulatedVector<Real> {
private:
  const std::vector<ROL::Ptr<Vector<Real> > >    vecs_;
  const ROL::Ptr<BatchManager<Real> >            bman_;
  const ROL::Ptr<SampleGenerator<Real> >         sampler_;
  mutable std::vector<ROL::Ptr<Vector<Real> > >  primal_vecs_;
  mutable ROL::Ptr<PrimalSimulatedVector<Real> > primal_pvec_;
  mutable bool isPrimalInitialized_;
public:

  DualSimulatedVector(const std::vector<ROL::Ptr<Vector<Real> > > &vecs,
                      const ROL::Ptr<BatchManager<Real> >         &bman,
                      const ROL::Ptr<SampleGenerator<Real> >      &sampler)
    : SimulatedVector<Real>(vecs,bman), vecs_(vecs), bman_(bman), sampler_(sampler),
      isPrimalInitialized_(false) {
    for( int i=0; i<sampler_->numMySamples(); ++i ) {
      primal_vecs_.push_back((vecs_[i]->dual()).clone());
    }
  }

  Real dot(const Vector<Real> &x) const {
    const SimulatedVector<Real> &xs
      = dynamic_cast<const SimulatedVector<Real>&>(x);

   ROL_TEST_FOR_EXCEPTION( sampler_->numMySamples() != static_cast<Real>(xs.numVectors()),
                               std::invalid_argument,
                               "Error: Vectors must have the same number of subvectors." );

    Real c = 0;
    Real locresult = 0;
    Real result = 0;
    for( int i=0; i<sampler_->numMySamples(); ++i ) {
      //locresult += vecs_[i]->dot(*xs.get(i)) / sampler_->getMyWeight(i);
      Real y = vecs_[i]->dot(*xs.get(i)) / sampler_->getMyWeight(i) - c;
      Real t = locresult + y;
      c = (t - locresult) - y;
      locresult = t;
    }

    bman_->sumAll(&locresult, &result, 1);

    return result;
  }

  ROL::Ptr<Vector<Real> > clone(void) const {
    std::vector<ROL::Ptr<Vector<Real> > > clonevec;
    for( int i=0; i<sampler_->numMySamples(); ++i ) {
      clonevec.push_back(vecs_[i]->clone());
    }
    return ROL::makePtr<DualSimulatedVector<Real>>(clonevec, bman_, sampler_);
  }

  const Vector<Real>& dual(void) const {
    if (!isPrimalInitialized_) {
      primal_pvec_ = ROL::makePtr<PrimalSimulatedVector<Real>>(primal_vecs_, bman_, sampler_);
      isPrimalInitialized_ = true;
    }
    const Real one(1);
    for( int i=0; i<sampler_->numMySamples(); ++i ) {
      primal_vecs_[i]->set(vecs_[i]->dual());
      primal_vecs_[i]->scale(one/sampler_->getMyWeight(i));
    }
    return *primal_pvec_;
  }

};

template<class Real>
ROL::Ptr<const Vector<Real> > 
CreateSimulatedVector( const ROL::Ptr<const Vector<Real> > &a, 
                       const ROL::Ptr<BatchManager<Real> > &bman ) {
  
  typedef ROL::Ptr<const Vector<Real> >      Vp;
  typedef const SimulatedVector<Real> PV;

  Vp temp[] = {a};
  return ROL::makePtr<PV>( std::vector<Vp>(temp, temp+1), bman );
}

template<class Real>
ROL::Ptr<Vector<Real> > 
CreateSimulatedVector( const ROL::Ptr<Vector<Real> > &a,
                       const ROL::Ptr<Vector<Real> > &b,
                       const ROL::Ptr<BatchManager<Real> > &bman ) {
  
  typedef ROL::Ptr<Vector<Real> >      Vp;
  typedef SimulatedVector<Real> PV;

  Vp temp[] = {a,b};
  return ROL::makePtr<PV>( std::vector<Vp>(temp, temp+2), bman );
}

template<class Real>
ROL::Ptr<const Vector<Real> > 
CreateSimulatedVector( const ROL::Ptr<const Vector<Real> > &a,
                       const ROL::Ptr<const Vector<Real> > &b,
                       const ROL::Ptr<BatchManager<Real> > &bman ) {
  
  
  typedef ROL::Ptr<const Vector<Real> >      Vp;
  typedef const SimulatedVector<Real> PV;

  Vp temp[] = {a,b};
  return ROL::makePtr<PV>( std::vector<Vp>(temp, temp+2), bman );
}

template<class Real>
ROL::Ptr<Vector<Real> > 
CreateSimulatedVector( const ROL::Ptr<Vector<Real> > &a,
                       const ROL::Ptr<Vector<Real> > &b,
                       const ROL::Ptr<Vector<Real> > &c,
                       const ROL::Ptr<BatchManager<Real> > &bman ) {
  
  
  typedef ROL::Ptr<Vector<Real> >      Vp;
  typedef SimulatedVector<Real> PV;

  Vp temp[] = {a,b,c};
  return ROL::makePtr<PV>( std::vector<Vp>(temp, temp+3), bman );
}

template<class Real>
ROL::Ptr<const Vector<Real> > 
CreateSimulatedVector( const ROL::Ptr<const Vector<Real> > &a,
                       const ROL::Ptr<const Vector<Real> > &b,
                       const ROL::Ptr<const Vector<Real> > &c,
                       const ROL::Ptr<BatchManager<Real> > &bman ) {
  
  
  typedef ROL::Ptr<const Vector<Real> >      Vp;
  typedef const SimulatedVector<Real> PV;

  Vp temp[] = {a,b,c};
  return ROL::makePtr<PV>( std::vector<Vp>(temp, temp+3), bman );
}

template<class Real>
ROL::Ptr<Vector<Real> > 
CreateSimulatedVector( const ROL::Ptr<Vector<Real> > &a,
                       const ROL::Ptr<Vector<Real> > &b,
                       const ROL::Ptr<Vector<Real> > &c,
                       const ROL::Ptr<Vector<Real> > &d,
                       const ROL::Ptr<BatchManager<Real> > &bman ) {
  
  
  typedef ROL::Ptr<Vector<Real> >      Vp;
  typedef SimulatedVector<Real> PV;

  Vp temp[] = {a,b,c,d};
  return ROL::makePtr<PV>( std::vector<Vp>(temp, temp+4), bman );
}

template<class Real>
ROL::Ptr<const Vector<Real> > 
CreateSimulatedVector( const ROL::Ptr<const Vector<Real> > &a,
                       const ROL::Ptr<const Vector<Real> > &b,
                       const ROL::Ptr<const Vector<Real> > &c,
                       const ROL::Ptr<const Vector<Real> > &d,
                       const ROL::Ptr<BatchManager<Real> > &bman ) {
  
  
  typedef ROL::Ptr<const Vector<Real> >      Vp;
  typedef const SimulatedVector<Real> PV;

  Vp temp[] = {a,b,c,d};
  return ROL::makePtr<PV>( std::vector<Vp>(temp, temp+4), bman );
}

} // namespace ROL

#endif // ROL_SIMULATED_VECTOR_H

