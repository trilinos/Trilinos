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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_PRODUCT_VECTOR_SPACE_HPP
#define THYRA_DEFAULT_PRODUCT_VECTOR_SPACE_HPP


#include "Thyra_DefaultProductVectorSpace_decl.hpp"
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_DefaultProductMultiVector.hpp"
#include "Thyra_ProductMultiVectorBase.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_dyn_cast.hpp"


namespace Thyra {


// Constructors/initializers/accessors


template<class Scalar>
DefaultProductVectorSpace<Scalar>::DefaultProductVectorSpace()
  : numBlocks_(-1), dim_(-1)
{}


template<class Scalar>
DefaultProductVectorSpace<Scalar>::DefaultProductVectorSpace(
  const ArrayView<const RCP<const VectorSpaceBase<Scalar> > > &vecSpaces_in
  )
  : numBlocks_(-1), dim_(-1)
{
  initialize(vecSpaces_in);
}


template<class Scalar>
void DefaultProductVectorSpace<Scalar>::initialize(
  const ArrayView<const RCP<const VectorSpaceBase<Scalar> > > &vecSpaces_in
  )
{

  //
  // Check preconditions and compute cached quantities
  //
  const int nBlocks = vecSpaces_in.size();
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( nBlocks == 0 );
#endif
  bool overallHasInCoreView = true;
  for (int k = 0; k < nBlocks; ++k) {
#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPTION(
      vecSpaces_in[k].get() == NULL, std::invalid_argument
      ,"Error, the smart pointer vecSpaces["<<k<<"] can not be NULL!"
      );
#endif
    if (!vecSpaces_in[k]->hasInCoreView()) overallHasInCoreView = false;
  }

  //
  // Setup private data members (should not throw an exception from here)
  //
  numBlocks_ = nBlocks;
  vecSpaces_ = Teuchos::rcp(new vecSpaces_t);
  *vecSpaces_ = vecSpaces_in;
  vecSpacesOffsets_ = Teuchos::rcp(new vecSpacesOffsets_t(nBlocks+1));
  (*vecSpacesOffsets_)[0] = 0;
  dim_ = 0;
  for( int k = 1; k <= nBlocks; ++k ) {
    const Ordinal dim_km1 = vecSpaces_in[k-1]->dim();
    (*vecSpacesOffsets_)[k] = (*vecSpacesOffsets_)[k-1] + dim_km1;
    dim_ += dim_km1;
  }
  isInCore_ = overallHasInCoreView;

}


template<class Scalar>
void DefaultProductVectorSpace<Scalar>::uninitialize(
  const ArrayView<RCP<const VectorSpaceBase<Scalar> > > &vecSpaces_in
  )
{
  TEST_FOR_EXCEPT(!is_null(vecSpaces_in)); // ToDo: Implement!
  vecSpaces_ = Teuchos::null;
  vecSpacesOffsets_ = Teuchos::null;
  numBlocks_ = -1;
  dim_ = -1;
  isInCore_ = false;
}


template<class Scalar>
void DefaultProductVectorSpace<Scalar>::getVecSpcPoss(
  Ordinal i, int* kth_vector_space, Ordinal* kth_global_offset
  ) const
{
  // Validate the preconditions
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(
    !(0 <= i && i < this->dim()), std::out_of_range
    ,"VectorSpaceBlocked::get_vector_space_position(...): Error, i = "
    << i << " is not in range [0,"<<(this->dim()-1)<<"]"
    );
#endif
  *kth_vector_space  = 0;
  *kth_global_offset = 0;
  while( *kth_vector_space < numBlocks_ ) {
    const Ordinal off_kp1 = (*vecSpacesOffsets_)[*kth_vector_space+1];
    if( off_kp1 > i ) {
      *kth_global_offset = (*vecSpacesOffsets_)[*kth_vector_space];
      break;
    }
    ++(*kth_vector_space);
  }
  TEST_FOR_EXCEPT( !(*kth_vector_space < numBlocks_) );
}


// Overridden from DefaultProductVectorSpace


template<class Scalar>
int DefaultProductVectorSpace<Scalar>::numBlocks() const
{
  return numBlocks_;
}


template<class Scalar>
Teuchos::RCP<const VectorSpaceBase<Scalar> >
DefaultProductVectorSpace<Scalar>::getBlock(const int k) const
{
  TEST_FOR_EXCEPT( k < 0 || numBlocks_ < k );
  return (*vecSpaces_)[k];
}


// Overridden from VectorSpaceBase


template<class Scalar>
Ordinal DefaultProductVectorSpace<Scalar>::dim() const
{
  return dim_;
}


template<class Scalar>
bool DefaultProductVectorSpace<Scalar>::isCompatible(
  const VectorSpaceBase<Scalar>& vecSpc ) const
{

  using Teuchos::ptrFromRef;
  using Teuchos::ptr_dynamic_cast;

  const int nBlocks = this->numBlocks();

  // Check for product vector interface
  const Ptr<const ProductVectorSpaceBase<Scalar> > pvsb =
    ptr_dynamic_cast<const ProductVectorSpaceBase<Scalar> >(ptrFromRef(vecSpc));

  if (nonnull(pvsb)) {
    // Validate that constituent vector spaces are compatible
    if( nBlocks != pvsb->numBlocks() )
      return false;
    for( int i = 0; i < nBlocks; ++i ) {
      if( !this->getBlock(i)->isCompatible(*pvsb->getBlock(i)) )
        return false;
    }
    return true;
  }

  // Check for a single vector single vector space
  if (nBlocks == 1) {
    return this->getBlock(0)->isCompatible(vecSpc);
  }

  // If we get here, the RHS is not a product vector space and/or this is not
  // a single block VS so we can assume the spaces are *not* compatible!
  return false;

}


template<class Scalar>
Teuchos::RCP< VectorBase<Scalar> >
DefaultProductVectorSpace<Scalar>::createMember() const
{
  return defaultProductVector<Scalar>(Teuchos::rcpFromRef(*this));
}


template<class Scalar>
Scalar DefaultProductVectorSpace<Scalar>::scalarProd(
  const VectorBase<Scalar> &x_in,
  const VectorBase<Scalar> &y_in
  ) const
{
  const int nBlocks = this->numBlocks(); 
  const ProductVectorBase<Scalar>
    &x = Teuchos::dyn_cast<const ProductVectorBase<Scalar> >(x_in),
    &y = Teuchos::dyn_cast<const ProductVectorBase<Scalar> >(y_in);
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(
    nBlocks!=x.productSpace()->numBlocks()
    || nBlocks!=y.productSpace()->numBlocks()
    );
#endif
  Scalar scalarProd_rtn = Teuchos::ScalarTraits<Scalar>::zero();
  for( int k = 0; k < nBlocks; ++k )
    scalarProd_rtn += (*vecSpaces_)[k]->scalarProd(
      *x.getVectorBlock(k),*y.getVectorBlock(k)
      );
  return scalarProd_rtn;
}


template<class Scalar>
void DefaultProductVectorSpace<Scalar>::scalarProdsImpl(
  const MultiVectorBase<Scalar> &X_in,
 const MultiVectorBase<Scalar> &Y_in,
  const ArrayView<Scalar> &scalarProds_out
  ) const
{
  using Teuchos::as;
  using Teuchos::Workspace;
  const VectorSpaceBase<Scalar> &domain = *X_in.domain();
  const Ordinal m = domain.dim();
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(is_null(scalarProds_out));
  TEST_FOR_EXCEPT( !domain.isCompatible(*Y_in.domain()) );
  TEUCHOS_ASSERT_EQUALITY( as<Ordinal>(scalarProds_out.size()),
    as<Ordinal>(m) )
#endif
  if(m==1) {
    scalarProds_out[0] = this->scalarProd(*X_in.col(0),*Y_in.col(0));
    return;
    // ToDo: Remove this if(...) block once we have a DefaultProductMultiVector implementation!
  }
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();
  const int nBlocks = this->numBlocks(); 
  const ProductMultiVectorBase<Scalar>
    &X = Teuchos::dyn_cast<const ProductMultiVectorBase<Scalar> >(X_in),
    &Y = Teuchos::dyn_cast<const ProductMultiVectorBase<Scalar> >(Y_in);
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( nBlocks!=X.productSpace()->numBlocks() || nBlocks!=Y.productSpace()->numBlocks() );
#endif
  Workspace<Scalar> _scalarProds_out(wss, m, false);
  std::fill( scalarProds_out.begin(), scalarProds_out.end(),
    ScalarTraits<Scalar>::zero() );
  for( int k = 0; k < nBlocks; ++k ) {
    (*vecSpaces_)[k]->scalarProds(
      *X.getMultiVectorBlock(k), *Y.getMultiVectorBlock(k), _scalarProds_out());
    for( int j = 0; j < m; ++j )
      scalarProds_out[j] += _scalarProds_out[j];
  }
}


template<class Scalar>
bool DefaultProductVectorSpace<Scalar>::hasInCoreView(const Range1D& rng_in, const EViewType viewType, const EStrideType strideType) const
{
  const Range1D rng = full_range(rng_in,0,dim_-1);
  // First see if rng fits in a single constituent vector
  int    kth_vector_space  = -1;
  Ordinal  kth_global_offset = 0;
  this->getVecSpcPoss(rng.lbound(),&kth_vector_space,&kth_global_offset);
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( !( 0 <= kth_vector_space && kth_vector_space <= numBlocks_ ) );
#endif
  if( rng.lbound() + rng.size() <= kth_global_offset + (*vecSpaces_)[kth_vector_space]->dim() ) {
    return (*vecSpaces_)[kth_vector_space]->hasInCoreView(rng_in-kth_global_offset,viewType,strideType);
  }
  // If we get here, rng does not fit in a single constituent vector which
  // also means that numBlocks_ > 1 must also be true!
  //
  // Next, if the client is asking for a direct view then we have to return
  // false since this range spans more than one constituent vector.
  if( viewType == VIEW_TYPE_DIRECT )
    return false;
  // If we get here then hasDirectView==false and therefore we are allowed to
  // create a copy.  Therefore, if all of the constituent vectors are "in
  // core" then we can return true.
  if(isInCore_)
    return true;
  // Finally, loop through all of the constituent vectors spaned by rng and
  // see if they are each in core.
  //
  // Todo: Implement this if you have to!
  //
  // We must give up and return false
  return false;
}


template<class Scalar>
Teuchos::RCP< const VectorSpaceFactoryBase<Scalar> >
DefaultProductVectorSpace<Scalar>::smallVecSpcFcty() const
{
  if (dim_)
    return (*vecSpaces_)[0]->smallVecSpcFcty(); // They should all be compatible?
  return Teuchos::null;
}


template<class Scalar>
Teuchos::RCP< MultiVectorBase<Scalar> >
DefaultProductVectorSpace<Scalar>::createMembers(int numMembers) const
{
  return defaultProductMultiVector<Scalar>(Teuchos::rcpFromRef(*this),
    numMembers);
}


template<class Scalar>
Teuchos::RCP< const VectorSpaceBase<Scalar> >
DefaultProductVectorSpace<Scalar>::clone() const
{
  // Warning! If the client uninitialized this object then changes the
  // constituent vector spaces then we are in trouble!  The client is warned
  // in documentation!
  Teuchos::RCP<DefaultProductVectorSpace<Scalar> >
    pvs = productVectorSpace<Scalar>();
  pvs->numBlocks_          = numBlocks_;
  pvs->vecSpaces_          = vecSpaces_;
  pvs->vecSpacesOffsets_   = vecSpacesOffsets_;
  pvs->dim_                = dim_;
  pvs->isInCore_           = isInCore_;
  return pvs;
}


// Overridden from Teuchos::Describable

                                                
template<class Scalar>
std::string DefaultProductVectorSpace<Scalar>::description() const
{
  std::ostringstream oss;
  oss
    << Teuchos::Describable::description() << "{"
    << "dim="<<dim_
    << ",numBlocks="<<numBlocks_
    << "}";
  return oss.str();
}


template<class Scalar>
void DefaultProductVectorSpace<Scalar>::describe(
  Teuchos::FancyOStream                &out_arg
  ,const Teuchos::EVerbosityLevel      verbLevel
  ) const
{
  typedef Teuchos::ScalarTraits<Scalar>  ST;
  using Teuchos::includesVerbLevel;
  using Teuchos::OSTab;
  RCP<FancyOStream> out = rcpFromRef(out_arg);
  OSTab tab(out);
  if (includesVerbLevel(verbLevel, Teuchos::VERB_LOW, true)) {
    *out << this->description() << std::endl;
  }
  if (includesVerbLevel(verbLevel, Teuchos::VERB_MEDIUM) && numBlocks_ > 0) {
    OSTab tab2(out);
    *out <<  "Constituent vector spaces V[0], V[1], ... V[numBlocks-1]:\n";
    OSTab tab3(out);
    for( int k = 0; k < numBlocks_; ++k ) {
      *out << "V["<<k<<"] = " << Teuchos::describe(*(*vecSpaces_)[k],verbLevel);
    }
  }
}


} // namespace Thyra


#endif // THYRA_DEFAULT_PRODUCT_VECTOR_SPACE_HPP
