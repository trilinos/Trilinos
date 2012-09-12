// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_VECTOR_DEFAULT_BASE_DEF_HPP
#define THYRA_VECTOR_DEFAULT_BASE_DEF_HPP


// Define to make some verbose output
//#define THYRA_VECTOR_VERBOSE_TO_ERROR_OUT


#include "Thyra_VectorDefaultBase_decl.hpp"
#include "Thyra_VectorSpaceFactoryBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorDefaultBase.hpp"
#include "Thyra_AssertOp.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "RTOpPack_ROpGetSubVector.hpp"
#include "RTOpPack_TOpSetSubVector.hpp"
#include "Teuchos_Assert.hpp"


#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
#  include "Teuchos_VerboseObject.hpp"
#  define THYRA_VECTOR_VERBOSE_OUT_STATEMENT \
     RCP<Teuchos::FancyOStream> dbgout = Teuchos::VerboseObjectBase::getDefaultOStream()
#endif // THYRA_VECTOR_VERBOSE_TO_ERROR_OUT



namespace Thyra {


// Overridden from Teuchos::Describable


template<class Scalar>
std::string VectorDefaultBase<Scalar>::description() const
{
  std::ostringstream oss;
  const RCP<const VectorSpaceBase<Scalar> > vs = this->space();
  oss << Teuchos::Describable::description();
  if(is_null(vs)) {
    oss << "{space=NULL}"; 
  }
  else {
    const Ordinal dim = vs->dim();
    oss << "{dim=" << dim << "}";
  }
  return oss.str();
}


template<class Scalar>
void VectorDefaultBase<Scalar>::describe(
  Teuchos::FancyOStream &out_arg,
  const Teuchos::EVerbosityLevel verbLevel
  ) const
{
  using Teuchos::FancyOStream;
  using Teuchos::OSTab;
  RCP<FancyOStream> out = Teuchos::rcpFromRef(out_arg);
  OSTab tab(out);
  *out << this->description() << "\n";
  if (this->space()->dim()) {
    tab.incrTab();
    if (verbLevel >= Teuchos::VERB_HIGH) {
      const ConstDetachedVectorView<Scalar> dvv(*this);
      for( Ordinal i = 0; i < dvv.subDim(); ++i )
        *out << i << ":" << dvv[i] << std::endl;
    }
  }
}


// Overridden from LinearOpBase


template<class Scalar>
RCP< const VectorSpaceBase<Scalar> >
VectorDefaultBase<Scalar>::range() const
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  THYRA_VECTOR_VERBOSE_OUT_STATEMENT;
  *dbgout << "\nThyra::VectorDefaultBase<"
          <<Teuchos::ScalarTraits<Scalar>::name()
          <<">::range() called!\n";
#endif
  return this->space();
}


template<class Scalar>
RCP< const VectorSpaceBase<Scalar> >
VectorDefaultBase<Scalar>::domain() const
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  THYRA_VECTOR_VERBOSE_OUT_STATEMENT;
  *dbgout << "\nThyra::VectorDefaultBase<"
          <<Teuchos::ScalarTraits<Scalar>::name()
          <<">::domain() called!\n";
#endif
  if(!domain_.get()) {
    domain_ = range()->smallVecSpcFcty()->createVecSpc(1);
  }
  return domain_;
}


// Overridden from MultiVectorBase


template<class Scalar>
RCP<MultiVectorBase<Scalar> >
VectorDefaultBase<Scalar>::clone_mv() const
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  THYRA_VECTOR_VERBOSE_OUT_STATEMENT;
  *dbgout << "\nThyra::VectorDefaultBase<"
          <<Teuchos::ScalarTraits<Scalar>::name()
          <<">::clone_mv() called!\n";
#endif
  return this->clone_v();
}


// Overridden from VectorBase


template<class Scalar>
RCP<VectorBase<Scalar> >
VectorDefaultBase<Scalar>::clone_v() const
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  THYRA_VECTOR_VERBOSE_OUT_STATEMENT;
  *dbgout << "\nThyra::VectorDefaultBase<"
          <<Teuchos::ScalarTraits<Scalar>::name()
          <<">::clone_v() called!\n";
#endif
  RCP<VectorBase<Scalar> > copy = createMember(this->space());
  assign(copy.ptr(), *this);
  return copy;
}


// protected


// Overridden protected functions from MultiVectorVectorBase


template<class Scalar>
RCP<VectorBase<Scalar> >
VectorDefaultBase<Scalar>::nonconstColImpl(Ordinal j)
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  THYRA_VECTOR_VERBOSE_OUT_STATEMENT;
  *dbgout << "\nThyra::VectorDefaultBase<"
          <<Teuchos::ScalarTraits<Scalar>::name()<<">::nonconstColImpl(j) called!\n";
#endif
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( j != 0 );
#endif
  return Teuchos::rcp(this,false);
}


template<class Scalar>
RCP<const MultiVectorBase<Scalar> >
VectorDefaultBase<Scalar>::contigSubViewImpl( const Range1D& col_rng ) const
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  THYRA_VECTOR_VERBOSE_OUT_STATEMENT;
  *dbgout << "\nThyra::VectorDefaultBase<"
          <<Teuchos::ScalarTraits<Scalar>::name()
          <<">::contigSubViewImpl(col_rng) const called!\n";
#endif
  validateColRng(col_rng);
  return Teuchos::rcp(this,false);
}


template<class Scalar>
RCP<MultiVectorBase<Scalar> >
VectorDefaultBase<Scalar>::nonconstContigSubViewImpl( const Range1D& col_rng )
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  THYRA_VECTOR_VERBOSE_OUT_STATEMENT;
  *dbgout << "\nThyra::VectorDefaultBase<"
          <<Teuchos::ScalarTraits<Scalar>::name()
          <<">::nonconstContigSubViewImpl(col_rng) called!\n";
#endif
  validateColRng(col_rng);
  return Teuchos::rcp(this,false);
}


template<class Scalar>
RCP<const MultiVectorBase<Scalar> >
VectorDefaultBase<Scalar>::nonContigSubViewImpl(
  const ArrayView<const int> &cols ) const
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  THYRA_VECTOR_VERBOSE_OUT_STATEMENT;
  *dbgout << "\nThyra::VectorDefaultBase<"
          <<Teuchos::ScalarTraits<Scalar>::name()
          <<">::nonContigSubViewImpl(cols) called!\n";
#endif
  validateColIndexes(cols);
  return Teuchos::rcp(this,false);
}


template<class Scalar>
RCP<MultiVectorBase<Scalar> >
VectorDefaultBase<Scalar>::nonconstNonContigSubViewImpl(
  const ArrayView<const int> &cols )
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  THYRA_VECTOR_VERBOSE_OUT_STATEMENT;
  *dbgout << "\nThyra::VectorDefaultBase<"
          <<Teuchos::ScalarTraits<Scalar>::name()
          <<">::nonconstNonContigSubViewImpl(cols) called!\n";
#endif
  validateColIndexes(cols);
  return Teuchos::rcp(this,false);
}


template<class Scalar>
void VectorDefaultBase<Scalar>::acquireDetachedMultiVectorViewImpl(
  const Range1D &rowRng,
  const Range1D &colRng,
  RTOpPack::ConstSubMultiVectorView<Scalar> *sub_mv
  ) const
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  THYRA_VECTOR_VERBOSE_OUT_STATEMENT;
  *dbgout << "\nThyra::VectorDefaultBase<"
          <<Teuchos::ScalarTraits<Scalar>::name()
          <<">::acquireDetachedMultiVectorViewImpl() const called!\n";
#endif
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(sub_mv==NULL);
#endif
  validateColRng(colRng);
  RTOpPack::ConstSubVectorView<Scalar> sv;
  this->acquireDetachedView(rowRng,&sv);
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( sv.stride() != 1 ); // Can't handle non-unit stride yet but we could
#endif
  sub_mv->initialize( sv.globalOffset(), sv.subDim(), 0, 1, sv.values(), sv.subDim() );
}


template<class Scalar>
void VectorDefaultBase<Scalar>::releaseDetachedMultiVectorViewImpl(
  RTOpPack::ConstSubMultiVectorView<Scalar>* sub_mv
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPT(sub_mv == 0);
  sub_mv->uninitialize();
}


template<class Scalar>
void VectorDefaultBase<Scalar>::acquireNonconstDetachedMultiVectorViewImpl(
  const Range1D &rowRng,
  const Range1D &colRng,
  RTOpPack::SubMultiVectorView<Scalar> *sub_mv
  )
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  THYRA_VECTOR_VERBOSE_OUT_STATEMENT;
  *dbgout << "\nThyra::VectorDefaultBase<"
          <<Teuchos::ScalarTraits<Scalar>::name()
          <<">::acquireNonconstDetachedMultiVectorViewImpl() called!\n";
#endif
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(sub_mv==NULL);
#endif
  validateColRng(colRng);
  RTOpPack::SubVectorView<Scalar> sv;
  this->acquireDetachedView(rowRng,&sv);
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( sv.stride() != 1 ); // Can't handle non-unit stride yet but we could
#endif
  sub_mv->initialize( sv.globalOffset(), sv.subDim(), 0, 1, sv.values(), sv.subDim() );
}


template<class Scalar>
void VectorDefaultBase<Scalar>::commitNonconstDetachedMultiVectorViewImpl(
  RTOpPack::SubMultiVectorView<Scalar>* sub_mv
  )
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  THYRA_VECTOR_VERBOSE_OUT_STATEMENT;
  *dbgout << "\nThyra::VectorDefaultBase<"
          <<Teuchos::ScalarTraits<Scalar>::name()
          <<">::commitNonconstDetachedMultiVectorViewImpl() called!\n";
#endif
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(sub_mv==NULL);
#endif
  RTOpPack::SubVectorView<Scalar> sv(
    sub_mv->globalOffset(),sub_mv->subDim(),sub_mv->values(),1);
  this->commitDetachedView(&sv);
  sub_mv->uninitialize();
}


// Overridden protected functions from VectorBase


template<class Scalar>
void VectorDefaultBase<Scalar>::acquireDetachedVectorViewImpl(
  const Range1D& rng_in, RTOpPack::ConstSubVectorView<Scalar>* sub_vec_inout
  ) const
{
  using Teuchos::dyn_cast;
  const Range1D rng = rng_in.full_range() ? Range1D(0,this->space()->dim()-1) : rng_in;
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(rng.ubound() < this->space()->dim()), std::out_of_range
    ,"VectorDefaultBase<Scalar>::acquireDetachedVectorViewImpl(rng,...):"
    " Error, rng = ["<<rng.lbound()<<","<<rng.ubound()
    <<"] is not in range = [0,"<<(this->space()->dim()-1)<<"]" );
#endif
  // Initialize the operator
  RTOpPack::ROpGetSubVector<Scalar> get_sub_vector_op(rng.lbound(),rng.ubound());
  // Create the reduction object (another sub_vec)
  RCP<RTOpPack::ReductTarget>
    reduct_obj = get_sub_vector_op.reduct_obj_create(); // This is really of type RTOpPack::ConstSubVectorView<Scalar>!
  // Perform the reduction (get the sub-vector requested)
  ::Thyra::applyOp<Scalar>(get_sub_vector_op, tuple(Teuchos::ptr<const VectorBase<Scalar> >(this))(),
    Teuchos::null, reduct_obj.ptr());
  // Get the sub-vector.
  *sub_vec_inout = get_sub_vector_op(*reduct_obj);
}


template<class Scalar>
void VectorDefaultBase<Scalar>::releaseDetachedVectorViewImpl(
  RTOpPack::ConstSubVectorView<Scalar>* sub_vec
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPT(sub_vec == 0);
  sub_vec->uninitialize();
}


template<class Scalar>
void VectorDefaultBase<Scalar>::acquireNonconstDetachedVectorViewImpl(
  const Range1D& rng, RTOpPack::SubVectorView<Scalar>* sub_vec_inout
  )
{
  //
  // Here we get a copy of the data for the sub-vector that the
  // client will modify.  We must later commit these changes to the
  // actual vector when the client calls commitDetachedView(...).
  // Note, this implementation is very dependent on the behavior of
  // the default implementation of constant version of
  // VectorDefaultBase<Scalar>::acquireDetachedView(...) and the implementation of
  // VectorDefaultBase<Scalar>::setSubVector(...)!
  //
  RTOpPack::ConstSubVectorView<Scalar> sub_vec;
  VectorDefaultBase<Scalar>::acquireDetachedVectorViewImpl( rng, &sub_vec );
  sub_vec_inout->initialize(
    sub_vec.globalOffset(), sub_vec.subDim(),
    Teuchos::arcp_const_cast<Scalar>(sub_vec.values()), sub_vec.stride()
    );
}


template<class Scalar>
void VectorDefaultBase<Scalar>::commitNonconstDetachedVectorViewImpl(
  RTOpPack::SubVectorView<Scalar>* sub_vec_inout
  )
{
  TEUCHOS_TEST_FOR_EXCEPT(sub_vec_inout == 0);
  RTOpPack::SparseSubVectorT<Scalar> spc_sub_vec(
    sub_vec_inout->globalOffset(), sub_vec_inout->subDim()
    ,sub_vec_inout->values(), sub_vec_inout->stride()
    );
  VectorDefaultBase<Scalar>::setSubVectorImpl(spc_sub_vec); // Commit the changes!
  sub_vec_inout->uninitialize(); // Make null as promised!
}


template<class Scalar>
void VectorDefaultBase<Scalar>::setSubVectorImpl( const RTOpPack::SparseSubVectorT<Scalar>& sub_vec )
{
  RTOpPack::TOpSetSubVector<Scalar> set_sub_vector_op(sub_vec);
  ::Thyra::applyOp<Scalar>(set_sub_vector_op, Teuchos::null,
    Teuchos::tuple(Teuchos::ptr<VectorBase<Scalar> >(this))(), Teuchos::null);
}


// Overridden protected functions from LinearOpBase


template<class Scalar>
bool VectorDefaultBase<Scalar>::opSupportedImpl(EOpTransp M_trans) const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  return ( ST::isComplex ? ( M_trans==NOTRANS || M_trans==CONJTRANS ) : true );
}


template<class Scalar>
void VectorDefaultBase<Scalar>::applyImpl(
  const EOpTransp M_trans,
  const MultiVectorBase<Scalar> &X,
  const Ptr<MultiVectorBase<Scalar> > &Y,
  const Scalar alpha,
  const Scalar beta
  ) const
{

  typedef Teuchos::ScalarTraits<Scalar> ST;

  // Validate input
#ifdef TEUCHOS_DEBUG
  THYRA_ASSERT_LINEAR_OP_MULTIVEC_APPLY_SPACES(
    "VectorDefaultBase<Scalar>::apply()", *this, M_trans, X, &*Y);
#endif

  const Ordinal numCols = X.domain()->dim();

  for (Ordinal col_j = 0; col_j < numCols; ++col_j) {

    // Get single column vectors
    const RCP<const VectorBase<Scalar> > x = X.col(col_j);
    const RCP<VectorBase<Scalar> > y = Y->col(col_j);

    // Here M = m (where m is a column vector)
    if( M_trans == NOTRANS || (M_trans == CONJ && !ST::isComplex) ) {
      // y = beta*y + alpha*m*x  (x is a scalar!)
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
      THYRA_VECTOR_VERBOSE_OUT_STATEMENT;
      *dbgout << "\nThyra::VectorDefaultBase<"
              <<Teuchos::ScalarTraits<Scalar>::name()
              <<">::apply(...) : y = beta*y + alpha*m*x  (x is a scalar!)\n";
#endif
      Vt_S( y.ptr(), beta );
      Vp_StV( y.ptr(), Scalar(alpha*get_ele(*x,0)), *this );
    }
    else if( M_trans == CONJTRANS || (M_trans == TRANS && !ST::isComplex) ) {
      // y = beta*y + alpha*m'*x  (y is a scalar!)
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
      THYRA_VECTOR_VERBOSE_OUT_STATEMENT;
      *dbgout << "\nThyra::VectorDefaultBase<"
              <<Teuchos::ScalarTraits<Scalar>::name()
              <<">::apply(...) : y = beta*y + alpha*m'*x  (y is a scalar!)\n";
#endif
      Scalar y_inout;
      if( beta == ST::zero() ) {
        y_inout = ST::zero();
      }
      else {
        y_inout = beta*get_ele(*y,0);
      }
#if defined(THYRA_VECTOR_VERBOSE_TO_ERROR_OUT) && defined(RTOPPACK_SPMD_APPLY_OP_DUMP)
      RTOpPack::show_spmd_apply_op_dump = true;
#endif
#if defined(THYRA_VECTOR_VERBOSE_TO_ERROR_OUT) && defined(RTOPPACK_RTOPT_HELPER_DUMP_OUTPUT)
      RTOpPack::rtop_helpers_dump_all = true;
#endif
      y_inout += alpha * this->space()->scalarProd(*this, *x);
#if defined(THYRA_VECTOR_VERBOSE_TO_ERROR_OUT) && defined(RTOPPACK_SPMD_APPLY_OP_DUMP)
      RTOpPack::show_spmd_apply_op_dump = false;
#endif
#if defined(THYRA_VECTOR_VERBOSE_TO_ERROR_OUT) && defined(RTOPPACK_RTOPT_HELPER_DUMP_OUTPUT)
      RTOpPack::rtop_helpers_dump_all = false;
#endif
      set_ele(0, y_inout, y.ptr());
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
      *dbgout
        << "\nThyra::VectorDefaultBase<"<<ST::name()<<">::apply(...) : y_inout = "
        << y_inout << "\n";
#endif
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "VectorBase<"<<ST::name()<<">::apply(M_trans,...): Error, M_trans="
        <<toString(M_trans)<<" not supported!" );
    }
    
  }

}


// private


template<class Scalar>
inline
void VectorDefaultBase<Scalar>::validateColRng( const Range1D &col_rng ) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(
    !( col_rng.full_range() || ( col_rng.lbound() == 0 && col_rng.ubound() == 0) ) );
#endif
}


template<class Scalar>
inline
void VectorDefaultBase<Scalar>::validateColIndexes(
  const ArrayView<const int>&cols ) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( cols.size() != 1 || cols[0] != 0 );
#endif
}


} // end namespace Thyra


#endif  // THYRA_VECTOR_DEFAULT_BASE_DEF_HPP
