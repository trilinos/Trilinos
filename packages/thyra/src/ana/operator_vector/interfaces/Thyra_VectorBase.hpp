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

#ifndef THYRA_VECTOR_HPP
#define THYRA_VECTOR_HPP

// Define to make some verbose output
//#define THYRA_VECTOR_VERBOSE_TO_ERROR_OUT

#include "Thyra_VectorBaseDecl.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_AssertOp.hpp"
#include "Thyra_SerialVectorSpaceStd.hpp"
#include "RTOpPack_ROpGetSubVector.hpp"
#include "RTOpPack_TOpSetSubVector.hpp"
#include "Teuchos_TestForException.hpp"

namespace Thyra {

template<class Scalar>
Teuchos::RefCountPtr<VectorBase<Scalar> >
VectorBase<Scalar>::clone_v() const
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nVector<Scalar>::clone_v() called!\n";
#endif
  Teuchos::RefCountPtr<VectorBase<Scalar> > copy = createMember(this->space());
  assign( &*copy, *this );
  return copy;
}

template<class Scalar>
void VectorBase<Scalar>::getSubVector( const Range1D& rng_in, RTOpPack::SubVectorT<Scalar>* sub_vec_inout ) const
{
  using Teuchos::dyn_cast;
  const Range1D rng = rng_in.full_range() ? Range1D(1,this->space()->dim()) : rng_in;
#ifdef _DEBUG
  TEST_FOR_EXCEPTION(
    this->space()->dim() < rng.ubound(), std::out_of_range
    ,"VectorBase<Scalar>::getSubVector(rng,...): Error, rng = ["<<rng.lbound()<<","<<rng.ubound()
    <<"] is not in range = [1,"<<this->space()->dim()<<"]" );
#endif
  // Free sub_vec if needed (note this is dependent on the implementation of this operator class!)
  if( sub_vec_inout->values() ) {
    free( (void*)sub_vec_inout->values()  );
  }
  // Initialize the operator
  RTOpPack::ROpGetSubVector<Scalar> get_sub_vector_op(rng.lbound(),rng.ubound());
  // Create the reduction object (another sub_vec)
  Teuchos::RefCountPtr<RTOpPack::ReductTarget>
    reduct_obj = get_sub_vector_op.reduct_obj_create(); // This is really of type RTOpPack::SubVectorT<Scalar>!
  // Perform the reduction (get the sub-vector requested)
  const VectorBase* sub_vecs[] = { this };
  applyOp(
    get_sub_vector_op,1,sub_vecs,0,NULL,&*reduct_obj
    ,rng.lbound(),rng.size(),rng.lbound()-1 // first_ele, sub_dim, global_offset
    );
  // Get the sub-vector.  Note reduct_obj will go out of scope so the sub_vec parameter will
  // own the memory allocated within get_sub_vector_op.create_reduct_obj_raw(...).  This is okay
  // since the client is required to call release_sub_vector(...) so release memory!
  RTOpPack::ReductTargetSubVectorT<Scalar> &sub_vec_ro = get_sub_vector_op(*reduct_obj);
  sub_vec_ro.transfer(sub_vec_inout);
}

template<class Scalar>
void VectorBase<Scalar>::freeSubVector( RTOpPack::SubVectorT<Scalar>* sub_vec ) const
{
  // Free sub_vec if needed (note this is dependent on the implementation of this operator class!)
  RTOpPack::ReductTargetSubVectorT<Scalar>::free(sub_vec);
}

template<class Scalar>
void VectorBase<Scalar>::getSubVector( const Range1D& rng, RTOpPack::MutableSubVectorT<Scalar>* sub_vec_inout )
{
  //
  // Here we get a copy of the data for the sub-vector that the
  // client will modify.  We must later commit these changes to the
  // actual vector when the client calls commitSubVector(...).
  // Note, this implementation is very dependent on the behavior of
  // the default implementation of constant version of
  // VectorBase<Scalar>::getSubVector(...) and the implementation of
  // VectorBase<Scalar>::setSubVector(...)!
  //
  RTOpPack::SubVectorT<Scalar> sub_vec;
  VectorBase<Scalar>::getSubVector( rng, &sub_vec );
  sub_vec_inout->initialize(
    sub_vec.globalOffset(),sub_vec.subDim(),const_cast<Scalar*>(sub_vec.values()),sub_vec.stride());
}

template<class Scalar>
void VectorBase<Scalar>::commitSubVector( RTOpPack::MutableSubVectorT<Scalar>* sub_vec_inout )
{
  RTOpPack::SparseSubVectorT<Scalar> spc_sub_vec(
    sub_vec_inout->globalOffset(), sub_vec_inout->subDim()
    ,sub_vec_inout->values(), sub_vec_inout->stride()
    );
  VectorBase<Scalar>::setSubVector( spc_sub_vec );       // Commit the changes!
  RTOpPack::SubVectorT<Scalar> sub_vec(*sub_vec_inout);
  VectorBase<Scalar>::freeSubVector( &sub_vec );         // Free the memory!
  sub_vec_inout->set_uninitialized();                    // Make null as promised!
}

template<class Scalar>
void VectorBase<Scalar>::setSubVector( const RTOpPack::SparseSubVectorT<Scalar>& sub_vec )
{
  RTOpPack::TOpSetSubVector<Scalar> set_sub_vector_op(sub_vec);
  VectorBase* targ_vecs[1] = { this };
  this->applyOp(
    set_sub_vector_op,0,NULL,1,targ_vecs,NULL
    ,sub_vec.globalOffset()+1,sub_vec.subDim(),sub_vec.globalOffset() // first_ele, sub_dim, global_offset
    );
}

// Overridden from Teuchos::Describable

template<class Scalar>
std::ostream& VectorBase<Scalar>::describe(
    std::ostream                         &out
    ,const Teuchos::EVerbosityLevel      verbLevel
    ,const std::string                   leadingIndent
    ,const std::string                   indentSpacer
    ) const
{
  out << leadingIndent << indentSpacer << "type = \'" << this->describe()
      << "\', size = " << this->space()->dim() << "\n";
  if(verbLevel >= Teuchos::VERB_HIGH) {
    RTOpPack::SubVectorT<Scalar> sv;
    this->getSubVector(Range1D(),&sv);
    for( Index i = 1; i <= sv.subDim(); ++i )
      out << leadingIndent << indentSpacer << indentSpacer << i << ":" << sv(i) << std::endl;
    this->freeSubVector(&sv);
  }
  return out;
}

// Overridden from OpBase

template<class Scalar>
bool VectorBase<Scalar>::opSupported(ETransp M_trans) const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  return ( ST::isComplex ? ( M_trans==NOTRANS || M_trans==CONJTRANS ) : true );
}

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> >
VectorBase<Scalar>::range() const
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nVector<Scalar>::range() called!\n";
#endif
  return this->space();
}

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> >
VectorBase<Scalar>::domain() const
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nVector<Scalar>::domain() called!\n";
#endif
  if(!domain_.get())
    const_cast<VectorBase<Scalar>*>(this)->domain_ = Teuchos::rcp(new SerialVectorSpaceStd<Scalar>(1));
  return domain_;
}

// Overridden from LinearOpBase

template<class Scalar>
void VectorBase<Scalar>::apply(
  const ETransp                M_trans
  ,const VectorBase<Scalar>    &x
  ,VectorBase<Scalar>          *y
  ,const Scalar                alpha
  ,const Scalar                beta
  ) const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  // Validate input
#ifdef _DEBUG
  THYRA_ASSERT_LINEAR_OP_VEC_APPLY_SPACES("VectorBase<Scalar>::apply()",*this,M_trans,x,y);
#endif
  // Here M = m (where m is a column vector)
  if( M_trans == NOTRANS || (M_trans == CONJ && !ST::isComplex) ) {
    // y = beta*y + alpha*m*x  (x is a scalar!)
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
    std::cerr << "\nVector<Scalar>::apply(...) : y = beta*y + alpha*m*x  (x is a scalar!)\n";
#endif
    Vt_S( y, beta );
    Vp_StV( y, Scalar(alpha*get_ele(x,1)), *this );
  }
  else if( M_trans == CONJTRANS || (M_trans == TRANS && !ST::isComplex) ) {
    // y = beta*y + alpha*m'*x  (y is a scalar!)
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
    std::cerr << "\nVector<Scalar>::apply(...) : y = beta*y + alpha*m'*x  (y is a scalar!)\n";
#endif
    Scalar y_inout;
    if( beta == ST::zero() ) {
      y_inout = ST::zero();
    }
    else {
      y_inout = get_ele(*y,1);
      y_inout = beta*y_inout;
    }
    y_inout += alpha*space()->scalarProd(*this,x);
    set_ele(1,y_inout,y);
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
void VectorBase<Scalar>::validateColRng( const Range1D &col_rng ) const
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT( !( col_rng.full_range() || ( col_rng.lbound() == 1 && col_rng.ubound() == 1) ) );
#endif
}

template<class Scalar>
inline
void VectorBase<Scalar>::validateColIndexes(  const int numCols, const int cols[] ) const
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT( numCols != 1 || cols == NULL || cols[0] != 1 );
#endif
}

template<class Scalar>
Teuchos::RefCountPtr<VectorBase<Scalar> >
VectorBase<Scalar>::col(Index j)
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nVector<Scalar>::col(j) called!\n";
#endif
#ifdef _DEBUG
  TEST_FOR_EXCEPT( j != 1 );
#endif
  return Teuchos::rcp(this,false);
}

template<class Scalar>
Teuchos::RefCountPtr<MultiVectorBase<Scalar> >
VectorBase<Scalar>::clone_mv() const
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nVector<Scalar>::clone_mv() called!\n";
#endif
  return this->clone_v();
}

template<class Scalar>
Teuchos::RefCountPtr<const MultiVectorBase<Scalar> >
VectorBase<Scalar>::subView( const Range1D& col_rng ) const
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nVector<Scalar>::subView(col_rng) const called!\n";
#endif
  validateColRng(col_rng);
  return Teuchos::rcp(this,false);
}

template<class Scalar>
Teuchos::RefCountPtr<MultiVectorBase<Scalar> >
VectorBase<Scalar>::subView( const Range1D& col_rng )
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nVector<Scalar>::subView(col_rng) called!\n";
#endif
  validateColRng(col_rng);
  return Teuchos::rcp(this,false);
}

template<class Scalar>
Teuchos::RefCountPtr<const MultiVectorBase<Scalar> >
VectorBase<Scalar>::subView( const int numCols, const int cols[] ) const
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nVector<Scalar>::subView(numCols,cols) const called!\n";
#endif
  validateColIndexes(numCols,cols);
  return Teuchos::rcp(this,false);
}

template<class Scalar>
Teuchos::RefCountPtr<MultiVectorBase<Scalar> >
VectorBase<Scalar>::subView( const int numCols, const int cols[] )
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nVector<Scalar>::subView(numCols,cols) called!\n";
#endif
  validateColIndexes(numCols,cols);
  return Teuchos::rcp(this,false);
}

template<class Scalar>
void VectorBase<Scalar>::getSubMultiVector(
  const Range1D                       &rowRng
  ,const Range1D                      &colRng
  ,RTOpPack::SubMultiVectorT<Scalar>  *sub_mv
  ) const
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nVector<Scalar>::getSubMultiVector() const called!\n";
#endif
#ifdef _DEBUG
  TEST_FOR_EXCEPT(sub_mv==NULL);
#endif
  validateColRng(colRng);
  RTOpPack::SubVectorT<Scalar> sv;
  getSubVector(rowRng,&sv);
#ifdef _DEBUG
  TEST_FOR_EXCEPT( sv.stride() != 1 ); // Can't handle non-unit stride yet but we could
#endif
  sub_mv->initialize( sv.globalOffset(), sv.subDim(), 0, 1, sv.values(), sv.subDim() );
}

template<class Scalar>
void VectorBase<Scalar>::freeSubMultiVector( RTOpPack::SubMultiVectorT<Scalar>* sub_mv ) const
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nVector<Scalar>::freeSubMultiVector() const called!\n";
#endif
#ifdef _DEBUG
  TEST_FOR_EXCEPT(sub_mv==NULL);
#endif
  RTOpPack::SubVectorT<Scalar> sv(sub_mv->globalOffset(),sub_mv->subDim(),sub_mv->values(),1);
  freeSubVector(&sv);
  sub_mv->set_uninitialized();
}

template<class Scalar>
void VectorBase<Scalar>::getSubMultiVector(
  const Range1D                                &rowRng
  ,const Range1D                               &colRng
  ,RTOpPack::MutableSubMultiVectorT<Scalar>    *sub_mv
  )
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nVector<Scalar>::getSubMultiVector() called!\n";
#endif
#ifdef _DEBUG
  TEST_FOR_EXCEPT(sub_mv==NULL);
#endif
  validateColRng(colRng);
  RTOpPack::MutableSubVectorT<Scalar> sv;
  getSubVector(rowRng,&sv);
#ifdef _DEBUG
  TEST_FOR_EXCEPT( sv.stride() != 1 ); // Can't handle non-unit stride yet but we could
#endif
  sub_mv->initialize( sv.globalOffset(), sv.subDim(), 0, 1, sv.values(), sv.subDim() );
}

template<class Scalar>
void VectorBase<Scalar>::commitSubMultiVector( RTOpPack::MutableSubMultiVectorT<Scalar>* sub_mv )
{
#ifdef THYRA_VECTOR_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nVector<Scalar>::commitSubMultiVector() called!\n";
#endif
#ifdef _DEBUG
  TEST_FOR_EXCEPT(sub_mv==NULL);
#endif
  RTOpPack::MutableSubVectorT<Scalar> sv(sub_mv->globalOffset(),sub_mv->subDim(),sub_mv->values(),1);
  commitSubVector(&sv);
  sub_mv->set_uninitialized();
}

} // end namespace Thyra

#endif  // THYRA_VECTOR_HPP
