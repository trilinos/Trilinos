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

#ifndef ROL_THYRAVECTOR_H
#define ROL_THYRAVECTOR_H

#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "ROL_Vector.hpp"

#include <exception>

/** \class ROL::ThyraVector
    \brief Implements the ROL::Vector interface for a Thyra Vector.
*/

namespace ROL {

template <class Real>
class ThyraVector : public Vector<Real> {
private:

  Teuchos::RCP<Thyra::VectorBase<Real> >  thyra_vec_;

  class GetEleAccelerator {
    Teuchos::RCP<const Thyra::VectorBase<Real> > vec_;

    std::vector<Teuchos::ArrayRCP<const Real> > flatVec_;

    GetEleAccelerator(); // hide default constructor

    // This function assumes Thyra vectors are just product and SPMD vectors and
    // produces a flat array of ArrayRCP objects
    std::vector<Teuchos::ArrayRCP<const Real> > buildFlatStructure(const Thyra::VectorBase<Real> & vec)
    { 
      using Teuchos::Ptr;
      using Teuchos::ptrFromRef;
      using Teuchos::ptr_dynamic_cast;

      // is an spmd vector?
      Ptr<const Thyra::SpmdVectorBase<Real> > spmd_vec = ptr_dynamic_cast<const Thyra::SpmdVectorBase<Real> >(ptrFromRef(vec));

      if(spmd_vec!=Teuchos::null) {
        std::vector<Teuchos::ArrayRCP<const Real> > flatVec(1);
        Teuchos::ArrayRCP<const Real> ptrArray;
        spmd_vec->getLocalData(ptrFromRef(ptrArray));

        flatVec[0] = ptrArray;

        return flatVec;
      }


      // it must be a product vector then
      Ptr<const Thyra::ProductVectorBase<Real> > prod_vec = ptr_dynamic_cast<const Thyra::ProductVectorBase<Real> >(ptrFromRef(vec));

      std::vector<Teuchos::ArrayRCP<const Real> > flatVec;
      for(int i=0;i<prod_vec->productSpace()->numBlocks();i++) {
        Teuchos::RCP<const Thyra::VectorBase<Real> > block = prod_vec->getVectorBlock(i);

        std::vector<Teuchos::ArrayRCP<const Real> > subVec = buildFlatStructure(*block);
        flatVec.insert(flatVec.end(),subVec.begin(),subVec.end());
      }

      return flatVec;
    }

  public: 
    GetEleAccelerator(const Teuchos::RCP<const Thyra::VectorBase<Real> > & vec)
      : vec_(vec) 
    { 
      flatVec_ = buildFlatStructure(*vec);
    }

    ::Thyra::Ordinal getSize() const { return vec_->space()->dim(); }
    ::Thyra::Ordinal getLocalSize() const 
    {
      Thyra::Ordinal sum=0;
      for(std::size_t b=0;b<flatVec_.size();b++)
        sum += flatVec_[b].size(); 
      return sum;
    }

    Real getEle(::Thyra::Ordinal i) const
    { return ::Thyra::get_ele(*vec_,i); }

    Real getLocalEle(::Thyra::Ordinal i) const
    { 
      Teuchos::Ordinal sum = 0;
      std::size_t b=0;
      for(b=0;b<flatVec_.size();b++) {
        if(i<sum+flatVec_[b].size())
          break;
        sum += flatVec_[b].size();
      }

      TEUCHOS_ASSERT(b<flatVec_.size());

      return flatVec_[b][i-sum];
    }
  };

  class SetGetEleAccelerator {
    Teuchos::RCP<Thyra::VectorBase<Real> > vec_;

    SetGetEleAccelerator(); // hide default constructor

    std::vector<Teuchos::ArrayRCP<Real> > flatVec_;

    // This function assumes Thyra vectors are just product and SPMD vectors and
    // produces a flat array of ArrayRCP objects
    std::vector<Teuchos::ArrayRCP<Real> > buildFlatStructure(Thyra::VectorBase<Real> & vec)
    { 
      using Teuchos::Ptr;
      using Teuchos::ptrFromRef;
      using Teuchos::ptr_dynamic_cast;

      // is an spmd vector?
      Ptr<Thyra::SpmdVectorBase<Real> > spmd_vec = ptr_dynamic_cast<Thyra::SpmdVectorBase<Real> >(ptrFromRef(vec));

      if(spmd_vec!=Teuchos::null) {
        std::vector<Teuchos::ArrayRCP<Real> > flatVec(1);
        Teuchos::ArrayRCP<Real> ptrArray;
        spmd_vec->getNonconstLocalData(Teuchos::ptrFromRef(ptrArray));

        flatVec[0] = ptrArray;
        return flatVec;
      }

      // it must be a product vector then
      Ptr<Thyra::ProductVectorBase<Real> > prod_vec = ptr_dynamic_cast<Thyra::ProductVectorBase<Real> >(ptrFromRef(vec));

      std::vector<Teuchos::ArrayRCP<Real> > flatVec;
      for(int i=0;i<prod_vec->productSpace()->numBlocks();i++) {
        Teuchos::RCP<Thyra::VectorBase<Real> > block = prod_vec->getNonconstVectorBlock(i);

        std::vector<Teuchos::ArrayRCP<Real> > subVec = buildFlatStructure(*block);
        flatVec.insert(flatVec.end(),subVec.begin(),subVec.end());
      }

      return flatVec;
    }
  public: 
    SetGetEleAccelerator(const Teuchos::RCP<Thyra::VectorBase<Real> > & vec)
      : vec_(vec) 
    { 
      flatVec_ = buildFlatStructure(*vec);
      flatVec_ = buildFlatStructure(*vec);
    }

    ::Thyra::Ordinal getSize() const { return vec_->space()->dim(); }
    ::Thyra::Ordinal getLocalSize() const 
    {
      Thyra::Ordinal sum=0;
      for(std::size_t b=0;b<flatVec_.size();b++)
        sum += flatVec_[b].size(); 
      return sum;
    }

    void setEle(::Thyra::Ordinal i,Real v) 
    { ::Thyra::set_ele(i,v,vec_.ptr()); }

    void setLocalEle(::Thyra::Ordinal i,Real v) const
    { 
      Teuchos::Ordinal sum = 0;
      std::size_t b=0;
      for(b=0;b<flatVec_.size();b++) {
        if(i<sum+flatVec_[b].size())
          break;
        sum += flatVec_[b].size();
      }

      TEUCHOS_ASSERT(b<flatVec_.size());

      flatVec_[b][i-sum] = v;
    }

    Real getEle(::Thyra::Ordinal i) const
    { return ::Thyra::get_ele(*vec_,i); }

    Real getLocalEle(::Thyra::Ordinal i) const
    { 
      Teuchos::Ordinal sum = 0;
      std::size_t b=0;
      for(b=0;b<flatVec_.size();b++) {
        if(i<sum+flatVec_[b].size())
          break;
        sum += flatVec_[b].size();
      }

     
      std::stringstream ss;
      ss << "Block identifier b= " << b << " is too large for i=" << i << " on array with " << getLocalSize() <<
            " and " << flatVec_.size() << " blocks.";
      ROL_TEST_FOR_EXCEPTION(b>=flatVec_.size(),std::logic_error, ss.str());

      return flatVec_[b][i-sum];
    }
  };

public:
  ~ThyraVector() {}

  ThyraVector(const Teuchos::RCP<Thyra::VectorBase<Real> > & thyra_vec) : thyra_vec_(thyra_vec) {}

  /** \brief Compute \f$y \leftarrow x + y\f$ where \f$y = \mbox{*this}\f$.
  */
  void plus( const Vector<Real> &x ) {
    const ThyraVector &ex = dynamic_cast<const ThyraVector&>(x);
    ::Thyra::Vp_V( thyra_vec_.ptr(), *ex.getVector());
  }

  /** \brief Compute \f$y \leftarrow \alpha y\f$ where \f$y = \mbox{*this}\f$.
  */
  void scale( const Real alpha ) { 
    ::Thyra::scale(alpha, thyra_vec_.ptr());
  }

  /** \brief Returns \f$ \langle y,x \rangle \f$ where \f$y = \mbox{*this}\f$.
  */
  Real dot( const Vector<Real> &x ) const {
    const ThyraVector &ex = dynamic_cast<const ThyraVector&>(x);
    return ::Thyra::dot<Real>(*thyra_vec_, *ex.thyra_vec_);
  }

  /** \brief Returns \f$ \| y \| \f$ where \f$y = \mbox{*this}\f$.
  */
  Real norm() const {
    // return ::Thyra::norm_2<Real>(*thyra_vec_);
    return std::sqrt(dot(*this));
  } 

  /** \brief Clone to make a new (uninitialized) vector.
  */
  Teuchos::RCP<Vector<Real> > clone() const{
    Teuchos::RCP<Thyra::VectorBase<Real> > tv = thyra_vec_->clone_v();
    return Teuchos::rcp( new ThyraVector(tv) ); 
  }

  /** \brief Compute \f$y \leftarrow \alpha x + y\f$ where \f$y = \mbox{*this}\f$.
  */
  void axpy( const Real alpha, const Vector<Real> &x ) {
    const ThyraVector &ex = dynamic_cast<const ThyraVector&>(x);
    ::Thyra::Vp_StV( thyra_vec_.ptr(), alpha, *ex.getVector());
  }

  /**  \brief Set to zero vector.
  */
  void zero() {
    ::Thyra::put_scalar(0.0, thyra_vec_.ptr());
  }

  /**  \brief Set all entries of the vector to alpha.
    */
  void setScalar(const Real C) {
      ::Thyra::put_scalar(C, thyra_vec_.ptr());
    }

  /**  \brief Set all entries of the vector to alpha.
    */
  void putScalar(Real alpha) {
      ::Thyra::put_scalar(alpha, thyra_vec_.ptr());
    }

  /**  \brief const get of the Thyra vector.
    */
  Teuchos::RCP<const Thyra::VectorBase<Real> > getVector() const {
    return thyra_vec_;
  }

  /**  \brief nonconst get of the Thyra vector.
    */
  Teuchos::RCP<Thyra::VectorBase<Real> > getVector()  {
    return thyra_vec_;
  }

  /** \brief Return i-th basis vector.
    */
  Teuchos::RCP<Vector<Real> > basis( const int i ) const {
    Teuchos::RCP<Vector<Real> > e = clone();
    Teuchos::RCP<Thyra::VectorBase<Real> > basisThyraVec = (Teuchos::rcp_static_cast<ThyraVector>(e))->getVector();
    ::Thyra::put_scalar(0.0, basisThyraVec.ptr());
    ::Thyra::set_ele(i,1.0, basisThyraVec.ptr());
    return e;
  }

  /** \brief Return dimension of the vector space.
    */
  int dimension() const {
    return thyra_vec_->space()->dim();
  }

  /** \brief Set \f$y \leftarrow x\f$ where \f$y = \mathtt{*this}\f$.
    */
  void set(const Vector<Real> &x ) {
    const ThyraVector &ex = dynamic_cast<const ThyraVector&>(x);
    ::Thyra::copy( *ex.getVector(), thyra_vec_.ptr() );
  }

  void applyUnary( const Elementwise::UnaryFunction<Real> &f ) {
    SetGetEleAccelerator thisAccel(thyra_vec_);

    for(::Thyra::Ordinal i=0;i<thisAccel.getLocalSize();i++) {
      Real val  = thisAccel.getLocalEle(i); 

      thisAccel.setLocalEle(i,f.apply(val));
    }
/*
    for(::Thyra::Ordinal i=0;i<thyra_vec_->space()->dim();i++) {
      Real val = ::Thyra::get_ele(*thyra_vec_,i);
      ::Thyra::set_ele(i,f.apply(val),thyra_vec_.ptr());
    }
*/
  }

  void applyBinary( const Elementwise::BinaryFunction<Real> &f, const Vector<Real> &x ) {
    const ThyraVector &ex = dynamic_cast<const ThyraVector&>(x);
    Teuchos::RCP< const Thyra::VectorBase<Real> > xp = ex.getVector();

    SetGetEleAccelerator thisAccel(thyra_vec_);
    GetEleAccelerator xpAccel(xp);

    TEUCHOS_ASSERT(thisAccel.getLocalSize()==xpAccel.getLocalSize());
    for(::Thyra::Ordinal i=0;i<thisAccel.getLocalSize();i++) {
      Real val  = thisAccel.getLocalEle(i); 
      Real xval = xpAccel.getLocalEle(i); 

      thisAccel.setLocalEle(i,f.apply(val,xval));
    }
/*
    TEUCHOS_ASSERT(thisAccel.getSize()==thyra_vec_->space()->dim());
    TEUCHOS_ASSERT(xpAccel.getSize()==xp->space()->dim());

    for(::Thyra::Ordinal i=0;i<thyra_vec_->space()->dim();i++) {
      Real val  = thisAccel.getEle(i); // ::Thyra::get_ele(*thyra_vec_,i);
      Real xval = xpAccel.getEle(i); // ::Thyra::get_ele(*xp,i);

      // ::Thyra::set_ele(i,f.apply(val,xval),thyra_vec_.ptr());
      thisAccel.setEle(i,f.apply(val,xval));
    }
*/
  }

  Real reduce( const Elementwise::ReductionOp<Real> &r ) const {
    Real result = r.initialValue();
    
    for(::Thyra::Ordinal i=0;i<thyra_vec_->space()->dim();i++) {
      r.reduce(::Thyra::get_ele(*thyra_vec_,i),result);  
    }
    return result;
  } 

}; // class ThyraVector

} // namespace ROL

#endif

