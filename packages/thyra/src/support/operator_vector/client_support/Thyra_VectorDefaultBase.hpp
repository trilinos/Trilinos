// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_VECTOR_DEFAULT_BASE_HPP
#define THYRA_VECTOR_DEFAULT_BASE_HPP

// Define to make some verbose output
//#define THYRA_VECTOR_VERBOSE_TO_ERROR_OUT

#include "Thyra_VectorDefaultBaseDecl.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_MultiVectorDefaultBase.hpp"
#include "Thyra_SingleRhsLinearOpBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_AssertOp.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "RTOpPack_ROpGetSubVector.hpp"
#include "RTOpPack_TOpSetSubVector.hpp"
#include "Teuchos_TestForException.hpp"

#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
#  include "Teuchos_VerboseObject.hpp"
#  define THYRA_VECTOR_VERBOSE_OUT_STATEMENT \
     Teuchos::RefCountPtr<Teuchos::FancyOStream> dbgout = Teuchos::VerboseObjectBase::getDefaultOStream()
#endif // THYRA_VECTOR_VERBOSE_TO_ERROR_OUT


namespace Thyra {


// Overridden from Teuchos::Describable


template<class Scalar>
std::string VectorDefaultBase<Scalar>::description() const
{
  std::ostringstream oss;
  const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >
    space = this->space();
  oss << Teuchos::Describable::description();
  if(is_null(space)) {
    oss << "{space=NULL}"; 
  }
  else {
    const Index dim = space->dim();
    oss << "{dim=" << dim << "}";
  }
  return oss.str();
}


template<class Scalar>
void VectorDefaultBase<Scalar>::describe(
    Teuchos::FancyOStream                &out_arg
    ,const Teuchos::EVerbosityLevel      verbLevel
    ) const
{
  using Teuchos::RefCountPtr;
  using Teuchos::FancyOStream;
  using Teuchos::OSTab;
  RefCountPtr<FancyOStream> out = rcp(&out_arg,false);
  OSTab tab(out);
  *out << this->description() << "\n";
  tab.incrTab();
  if (verbLevel >= Teuchos::VERB_HIGH) {
    RTOpPack::ConstSubVectorView<Scalar> sv;
    this->acquireDetachedView(Range1D(),&sv);
    for( Index i = 0; i < sv.subDim(); ++i )
      *out << i << ":" << sv(i) << std::endl;
    this->releaseDetachedView(&sv);
  }
}


// Overridden from LinearOpBase

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> >
VectorDefaultBase<Scalar>::range() const
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  THYRA_VECTOR_VERBOSE_OUT_STATEMENT;
  *dbgout << "\nThyra::VectorDefaultBase<"<<Teuchos::ScalarTraits<Scalar>::name()<<">::range() called!\n";
#endif
  return this->space();
}

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> >
VectorDefaultBase<Scalar>::domain() const
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  THYRA_VECTOR_VERBOSE_OUT_STATEMENT;
  *dbgout << "\nThyra::VectorDefaultBase<"<<Teuchos::ScalarTraits<Scalar>::name()<<">::domain() called!\n";
#endif
  if(!domain_.get())
    const_cast<VectorDefaultBase<Scalar>*>(this)->domain_
      = Teuchos::rcp(new DefaultSpmdVectorSpace<Scalar>(1));
  return domain_;
}

// Overridden from SingleRhsLinearOpBase

template<class Scalar>
bool VectorDefaultBase<Scalar>::opSupported(ETransp M_trans) const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  return ( ST::isComplex ? ( M_trans==NOTRANS || M_trans==CONJTRANS ) : true );
}

template<class Scalar>
void VectorDefaultBase<Scalar>::apply(
  const ETransp                M_trans
  ,const VectorBase<Scalar>    &x
  ,VectorBase<Scalar>          *y
  ,const Scalar                alpha
  ,const Scalar                beta
  ) const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  // Validate input
#ifdef TEUCHOS_DEBUG
  THYRA_ASSERT_LINEAR_OP_VEC_APPLY_SPACES("VectorDefaultBase<Scalar>::apply()",*this,M_trans,x,y);
#endif
  // Here M = m (where m is a column vector)
  if( M_trans == NOTRANS || (M_trans == CONJ && !ST::isComplex) ) {
    // y = beta*y + alpha*m*x  (x is a scalar!)
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
    THYRA_VECTOR_VERBOSE_OUT_STATEMENT;
    *dbgout << "\nThyra::VectorDefaultBase<"<<Teuchos::ScalarTraits<Scalar>::name()<<">::apply(...) : y = beta*y + alpha*m*x  (x is a scalar!)\n";
#endif
    Vt_S( y, beta );
    Vp_StV( y, Scalar(alpha*get_ele(x,0)), *this );
  }
  else if( M_trans == CONJTRANS || (M_trans == TRANS && !ST::isComplex) ) {
    // y = beta*y + alpha*m'*x  (y is a scalar!)
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
    THYRA_VECTOR_VERBOSE_OUT_STATEMENT;
    *dbgout << "\nThyra::VectorDefaultBase<"<<Teuchos::ScalarTraits<Scalar>::name()<<">::apply(...) : y = beta*y + alpha*m'*x  (y is a scalar!)\n";
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
    y_inout += alpha*this->space()->scalarProd(*this,x);
#if defined(THYRA_VECTOR_VERBOSE_TO_ERROR_OUT) && defined(RTOPPACK_SPMD_APPLY_OP_DUMP)
    RTOpPack::show_spmd_apply_op_dump = false;
#endif
#if defined(THYRA_VECTOR_VERBOSE_TO_ERROR_OUT) && defined(RTOPPACK_RTOPT_HELPER_DUMP_OUTPUT)
    RTOpPack::rtop_helpers_dump_all = false;
#endif
    set_ele(0,y_inout,y);
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
    *dbgout
      << "\nThyra::VectorDefaultBase<"<<ST::name()<<">::apply(...) : y_inout = "
      << y_inout << "\n";
#endif
  }
  else {
    TEST_FOR_EXCEPTION(
      true,std::logic_error
      ,"VectorBase<"<<ST::name()<<">::apply(M_trans,...): Error, M_trans="<<toString(M_trans)<<" not supported!"
      );
  }
}

// Overridden from MultiVectorBase

template<class Scalar>
inline
void VectorDefaultBase<Scalar>::validateColRng( const Range1D &col_rng ) const
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( !( col_rng.full_range() || ( col_rng.lbound() == 0 && col_rng.ubound() == 0) ) );
#endif
}

template<class Scalar>
inline
void VectorDefaultBase<Scalar>::validateColIndexes(  const int numCols, const int cols[] ) const
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( numCols != 1 || cols == NULL || cols[0] != 0 );
#endif
}

template<class Scalar>
Teuchos::RefCountPtr<VectorBase<Scalar> >
VectorDefaultBase<Scalar>::col(Index j)
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  THYRA_VECTOR_VERBOSE_OUT_STATEMENT;
  *dbgout << "\nThyra::VectorDefaultBase<"<<Teuchos::ScalarTraits<Scalar>::name()<<">::col(j) called!\n";
#endif
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( j != 0 );
#endif
  return Teuchos::rcp(this,false);
}

template<class Scalar>
Teuchos::RefCountPtr<MultiVectorBase<Scalar> >
VectorDefaultBase<Scalar>::clone_mv() const
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  THYRA_VECTOR_VERBOSE_OUT_STATEMENT;
  *dbgout << "\nThyra::VectorDefaultBase<"<<Teuchos::ScalarTraits<Scalar>::name()<<">::clone_mv() called!\n";
#endif
  return this->clone_v();
}

template<class Scalar>
Teuchos::RefCountPtr<const MultiVectorBase<Scalar> >
VectorDefaultBase<Scalar>::subView( const Range1D& col_rng ) const
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  THYRA_VECTOR_VERBOSE_OUT_STATEMENT;
  *dbgout << "\nThyra::VectorDefaultBase<"<<Teuchos::ScalarTraits<Scalar>::name()<<">::subView(col_rng) const called!\n";
#endif
  validateColRng(col_rng);
  return Teuchos::rcp(this,false);
}

template<class Scalar>
Teuchos::RefCountPtr<MultiVectorBase<Scalar> >
VectorDefaultBase<Scalar>::subView( const Range1D& col_rng )
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  THYRA_VECTOR_VERBOSE_OUT_STATEMENT;
  *dbgout << "\nThyra::VectorDefaultBase<"<<Teuchos::ScalarTraits<Scalar>::name()<<">::subView(col_rng) called!\n";
#endif
  validateColRng(col_rng);
  return Teuchos::rcp(this,false);
}

template<class Scalar>
Teuchos::RefCountPtr<const MultiVectorBase<Scalar> >
VectorDefaultBase<Scalar>::subView( const int numCols, const int cols[] ) const
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  THYRA_VECTOR_VERBOSE_OUT_STATEMENT;
  *dbgout << "\nThyra::VectorDefaultBase<"<<Teuchos::ScalarTraits<Scalar>::name()<<">::subView(numCols,cols) const called!\n";
#endif
  validateColIndexes(numCols,cols);
  return Teuchos::rcp(this,false);
}

template<class Scalar>
Teuchos::RefCountPtr<MultiVectorBase<Scalar> >
VectorDefaultBase<Scalar>::subView( const int numCols, const int cols[] )
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  THYRA_VECTOR_VERBOSE_OUT_STATEMENT;
  *dbgout << "\nThyra::VectorDefaultBase<"<<Teuchos::ScalarTraits<Scalar>::name()<<">::subView(numCols,cols) called!\n";
#endif
  validateColIndexes(numCols,cols);
  return Teuchos::rcp(this,false);
}

template<class Scalar>
void VectorDefaultBase<Scalar>::acquireDetachedView(
  const Range1D                       &rowRng
  ,const Range1D                      &colRng
  ,RTOpPack::ConstSubMultiVectorView<Scalar>  *sub_mv
  ) const
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  THYRA_VECTOR_VERBOSE_OUT_STATEMENT;
  *dbgout << "\nThyra::VectorDefaultBase<"<<Teuchos::ScalarTraits<Scalar>::name()<<">::acquireDetachedView() const called!\n";
#endif
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(sub_mv==NULL);
#endif
  validateColRng(colRng);
  RTOpPack::ConstSubVectorView<Scalar> sv;
  acquireDetachedView(rowRng,&sv);
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( sv.stride() != 1 ); // Can't handle non-unit stride yet but we could
#endif
  sub_mv->initialize( sv.globalOffset(), sv.subDim(), 0, 1, sv.values(), sv.subDim() );
}

template<class Scalar>
void VectorDefaultBase<Scalar>::releaseDetachedView( RTOpPack::ConstSubMultiVectorView<Scalar>* sub_mv ) const
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  THYRA_VECTOR_VERBOSE_OUT_STATEMENT;
  *dbgout << "\nThyra::VectorDefaultBase<"<<Teuchos::ScalarTraits<Scalar>::name()<<">::releaseDetachedView() const called!\n";
#endif
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(sub_mv==NULL);
#endif
  RTOpPack::ConstSubVectorView<Scalar> sv(sub_mv->globalOffset(),sub_mv->subDim(),sub_mv->values(),1);
  releaseDetachedView(&sv);
  sub_mv->set_uninitialized();
}

template<class Scalar>
void VectorDefaultBase<Scalar>::acquireDetachedView(
  const Range1D                                &rowRng
  ,const Range1D                               &colRng
  ,RTOpPack::SubMultiVectorView<Scalar>    *sub_mv
  )
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  THYRA_VECTOR_VERBOSE_OUT_STATEMENT;
  *dbgout << "\nThyra::VectorDefaultBase<"<<Teuchos::ScalarTraits<Scalar>::name()<<">::acquireDetachedView() called!\n";
#endif
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(sub_mv==NULL);
#endif
  validateColRng(colRng);
  RTOpPack::SubVectorView<Scalar> sv;
  acquireDetachedView(rowRng,&sv);
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( sv.stride() != 1 ); // Can't handle non-unit stride yet but we could
#endif
  sub_mv->initialize( sv.globalOffset(), sv.subDim(), 0, 1, sv.values(), sv.subDim() );
}

template<class Scalar>
void VectorDefaultBase<Scalar>::commitDetachedView( RTOpPack::SubMultiVectorView<Scalar>* sub_mv )
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  THYRA_VECTOR_VERBOSE_OUT_STATEMENT;
  *dbgout << "\nThyra::VectorDefaultBase<"<<Teuchos::ScalarTraits<Scalar>::name()<<">::commitDetachedView() called!\n";
#endif
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(sub_mv==NULL);
#endif
  RTOpPack::SubVectorView<Scalar> sv(sub_mv->globalOffset(),sub_mv->subDim(),sub_mv->values(),1);
  commitDetachedView(&sv);
  sub_mv->set_uninitialized();
}

// Overridden from VectorBase

template<class Scalar>
Teuchos::RefCountPtr<VectorBase<Scalar> >
VectorDefaultBase<Scalar>::clone_v() const
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  THYRA_VECTOR_VERBOSE_OUT_STATEMENT;
  *dbgout << "\nThyra::VectorDefaultBase<"<<Teuchos::ScalarTraits<Scalar>::name()<<">::clone_v() called!\n";
#endif
  Teuchos::RefCountPtr<VectorBase<Scalar> > copy = createMember(this->space());
  assign( &*copy, *this );
  return copy;
}

template<class Scalar>
void VectorDefaultBase<Scalar>::acquireDetachedView( const Range1D& rng_in, RTOpPack::ConstSubVectorView<Scalar>* sub_vec_inout ) const
{
  using Teuchos::dyn_cast;
  const Range1D rng = rng_in.full_range() ? Range1D(0,this->space()->dim()-1) : rng_in;
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(
    !(rng.ubound() < this->space()->dim()), std::out_of_range
    ,"VectorDefaultBase<Scalar>::acquireDetachedView(rng,...): Error, rng = ["<<rng.lbound()<<","<<rng.ubound()
    <<"] is not in range = [0,"<<(this->space()->dim()-1)<<"]" );
#endif
  // Free sub_vec if needed (note this is dependent on the implementation of this operator class!)
  if( sub_vec_inout->values() ) {
    free( (void*)sub_vec_inout->values()  );
  }
  // Initialize the operator
  RTOpPack::ROpGetSubVector<Scalar> get_sub_vector_op(rng.lbound(),rng.ubound());
  // Create the reduction object (another sub_vec)
  Teuchos::RefCountPtr<RTOpPack::ReductTarget>
    reduct_obj = get_sub_vector_op.reduct_obj_create(); // This is really of type RTOpPack::ConstSubVectorView<Scalar>!
  // Perform the reduction (get the sub-vector requested)
  const VectorBase<Scalar>* sub_vecs[] = { this };
  ::Thyra::applyOp<Scalar>(
    get_sub_vector_op,1,sub_vecs,0,NULL,&*reduct_obj
    ,rng.lbound(),rng.size(),rng.lbound() // first_ele_offset,sub_dim,global_offset
    );
  // Get the sub-vector.  Note reduct_obj will go out of scope so the sub_vec parameter will
  // own the memory allocated within get_sub_vector_op.create_reduct_obj_raw(...).  This is okay
  // since the client is required to call release_sub_vector(...) so release memory!
  RTOpPack::ReductTargetSubVectorT<Scalar> &sub_vec_ro = get_sub_vector_op(*reduct_obj);
  sub_vec_ro.transfer(sub_vec_inout);
}

template<class Scalar>
void VectorDefaultBase<Scalar>::releaseDetachedView( RTOpPack::ConstSubVectorView<Scalar>* sub_vec ) const
{
  // Free sub_vec if needed (note this is dependent on the implementation of this operator class!)
  RTOpPack::ReductTargetSubVectorT<Scalar>::free(sub_vec);
}

template<class Scalar>
void VectorDefaultBase<Scalar>::acquireDetachedView( const Range1D& rng, RTOpPack::SubVectorView<Scalar>* sub_vec_inout )
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
  VectorDefaultBase<Scalar>::acquireDetachedView( rng, &sub_vec );
  sub_vec_inout->initialize(
    sub_vec.globalOffset(),sub_vec.subDim(),const_cast<Scalar*>(sub_vec.values()),sub_vec.stride());
}

template<class Scalar>
void VectorDefaultBase<Scalar>::commitDetachedView( RTOpPack::SubVectorView<Scalar>* sub_vec_inout )
{
  RTOpPack::SparseSubVectorT<Scalar> spc_sub_vec(
    sub_vec_inout->globalOffset(), sub_vec_inout->subDim()
    ,sub_vec_inout->values(), sub_vec_inout->stride()
    );
  VectorDefaultBase<Scalar>::setSubVector( spc_sub_vec );       // Commit the changes!
  RTOpPack::ConstSubVectorView<Scalar> sub_vec(*sub_vec_inout);
  VectorDefaultBase<Scalar>::releaseDetachedView( &sub_vec );         // Free the memory!
  sub_vec_inout->set_uninitialized();                    // Make null as promised!
}

template<class Scalar>
void VectorDefaultBase<Scalar>::setSubVector( const RTOpPack::SparseSubVectorT<Scalar>& sub_vec )
{
  RTOpPack::TOpSetSubVector<Scalar> set_sub_vector_op(sub_vec);
  VectorBase<Scalar>* targ_vecs[1] = { this };
  ::Thyra::applyOp<Scalar>(
    set_sub_vector_op,0,NULL,1,targ_vecs,NULL
    ,sub_vec.globalOffset(),sub_vec.subDim(),sub_vec.globalOffset() // first_ele_offset,sub_dim,global_offset
    );
}

} // end namespace Thyra

#endif  // THYRA_VECTOR_DEFAULT_BASE_HPP
