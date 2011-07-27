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

#ifndef THYRA_DEFAULT_SPMD_MULTI_VECTOR_DEF_HPP
#define THYRA_DEFAULT_SPMD_MULTI_VECTOR_DEF_HPP

// Define to make some verbose output
//#define THYRA_DEFAULT_SPMD_MULTI_VECTOR_VERBOSE_TO_ERROR_OUT

#include "Thyra_DefaultSpmdMultiVector_decl.hpp"
#include "Thyra_SpmdMultiVectorBase.hpp"
#include "Thyra_VectorSpaceFactoryBase.hpp"
#include "Thyra_DefaultSpmdVector.hpp"
#include "Teuchos_Assert.hpp"


namespace Thyra {


//
// Simple utility class to copy back multi-vector entries
//
// ToDo: Refactor above to not use raw pointers!
//
template<class Scalar>
class CopyBackSpmdMultiVectorEntries {
public:
  CopyBackSpmdMultiVectorEntries(
    const ArrayView<const int> &cols,
    const ArrayRCP<const Scalar> &localValuesView, const Ordinal localSubDim,
    const ArrayRCP<Scalar> &localValues, const Ordinal leadingDim
    )
    : cols_(cols), localValuesView_(localValuesView), localSubDim_(localSubDim),
      localValues_(localValues), leadingDim_(leadingDim)
    {}
  ~CopyBackSpmdMultiVectorEntries()
    {
      typedef typename ArrayRCP<const Scalar>::const_iterator const_itr_t;
      typedef typename ArrayRCP<Scalar>::iterator itr_t;
      // Copy from contiguous storage column by column
      if (localValues_.strong_count()) {
        const int numCols = cols_.size();
        const const_itr_t lvv = localValuesView_.begin();
        const itr_t lv = localValues_.begin();
        for (int k = 0; k < numCols; ++k) {
          const int col_k = cols_[k];
          const const_itr_t lvv_k = lvv + localSubDim_*k;
          const itr_t lv_k = lv + leadingDim_*col_k;
          std::copy( lvv_k, lvv_k + localSubDim_, lv_k );
        }
      }
#ifdef THYRA_DEBUG
      else {
        ++DefaultSpmdMultiVector<Scalar>::numSkipCopyBack;
      }
#endif // THYRA_DEBUG
    }
private:
  Array<int> cols_;
  ArrayRCP<const Scalar> localValuesView_;
  Ordinal localSubDim_;
  ArrayRCP<Scalar> localValues_;
  Ordinal leadingDim_;
  // Not defined and not to be called
  CopyBackSpmdMultiVectorEntries();
  CopyBackSpmdMultiVectorEntries(const CopyBackSpmdMultiVectorEntries&);
  CopyBackSpmdMultiVectorEntries& operator=(const CopyBackSpmdMultiVectorEntries&);
};


template<class Scalar>
RCP<CopyBackSpmdMultiVectorEntries<Scalar> >
copyBackSpmdMultiVectorEntries(
  const ArrayView<const int> &cols,
  const ArrayRCP<const Scalar> &localValuesView, const Ordinal localSubDim,
  const ArrayRCP<Scalar> &localValues, const Ordinal leadingDim
  )
{
  return Teuchos::rcp(
    new CopyBackSpmdMultiVectorEntries<Scalar>(
      cols, localValuesView, localSubDim, localValues, leadingDim
      )
    );
}


//
// DefaultSpmdMultiVector
//


#ifdef THYRA_DEBUG
template<class Scalar>
int DefaultSpmdMultiVector<Scalar>::numSkipCopyBack(0);
#endif


// Constructors/initializers/accessors


template<class Scalar>
DefaultSpmdMultiVector<Scalar>::DefaultSpmdMultiVector()
  :leadingDim_(0)
{}


template<class Scalar>
DefaultSpmdMultiVector<Scalar>::DefaultSpmdMultiVector(
  const RCP<const SpmdVectorSpaceBase<Scalar> > &spmdRangeSpace
  ,const RCP<const ScalarProdVectorSpaceBase<Scalar> > &domainSpace
  )
{
  initialize(spmdRangeSpace,domainSpace);
}


template<class Scalar>
DefaultSpmdMultiVector<Scalar>::DefaultSpmdMultiVector(
  const RCP<const SpmdVectorSpaceBase<Scalar> > &spmdRangeSpace
  ,const RCP<const ScalarProdVectorSpaceBase<Scalar> > &domainSpace
  ,const ArrayRCP<Scalar> &localValues
  ,const Ordinal leadingDim
  )
{
  initialize(spmdRangeSpace,domainSpace,localValues,leadingDim);
}


template<class Scalar>
void DefaultSpmdMultiVector<Scalar>::initialize(
  const RCP<const SpmdVectorSpaceBase<Scalar> > &spmdRangeSpace
  ,const RCP<const ScalarProdVectorSpaceBase<Scalar> > &domainSpace
  )
{
  const Ordinal localSubDim = spmdRangeSpace->localSubDim();
  ArrayRCP<Scalar> values;
  if (localSubDim)
    values = Teuchos::arcp<Scalar>(localSubDim * domainSpace->dim());
  initialize(spmdRangeSpace, domainSpace, values, localSubDim);
}


template<class Scalar>
void DefaultSpmdMultiVector<Scalar>::initialize(
  const RCP<const SpmdVectorSpaceBase<Scalar> > &spmdRangeSpace,
  const RCP<const ScalarProdVectorSpaceBase<Scalar> > &domainSpace,
  const ArrayRCP<Scalar> &localValues,
  const Ordinal leadingDim_in
  )
{
  const Ordinal leadingDim =
    (leadingDim_in >= 0 ? leadingDim_in : spmdRangeSpace->localSubDim());
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT(!is_null(spmdRangeSpace));
  TEUCHOS_ASSERT(!is_null(domainSpace));
  if (spmdRangeSpace->dim()) {
    TEUCHOS_ASSERT(!is_null(localValues));
  }
  TEUCHOS_ASSERT_INEQUALITY(leadingDim, >=, spmdRangeSpace->localSubDim());
#endif
  spmdRangeSpace_ = spmdRangeSpace;
  domainSpace_ = domainSpace;
  localValues_ = localValues;
  leadingDim_ = leadingDim;
  this->updateSpmdSpace();
}


template<class Scalar>
void DefaultSpmdMultiVector<Scalar>::uninitialize(
  RCP<const SpmdVectorSpaceBase<Scalar> > *spmdRangeSpace
  ,RCP<const ScalarProdVectorSpaceBase<Scalar> > *domainSpace
  ,ArrayRCP<Scalar> *localValues
  ,Ordinal *leadingDim
  )
{
  if(spmdRangeSpace) *spmdRangeSpace = spmdRangeSpace_;
  if(domainSpace) *domainSpace = domainSpace_;
  if(localValues) *localValues = localValues_;
  if(leadingDim) *leadingDim = leadingDim_;

  spmdRangeSpace_ = Teuchos::null;
  domainSpace_ = Teuchos::null;
  localValues_ = Teuchos::null;
  leadingDim_ = 0;

  this->updateSpmdSpace();
}


template<class Scalar>
RCP< const ScalarProdVectorSpaceBase<Scalar> >
DefaultSpmdMultiVector<Scalar>::domainScalarProdVecSpc() const
{
#ifdef THYRA_DEFAULT_SPMD_MULTI_VECTOR_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nSpmdMultiVectorStd<Scalar>::domainScalarProdVecSpc() const called!\n";
#endif
  return domainSpace_;
}


// Overridden public functions from SpmdMultiVectorBase


template<class Scalar>
RCP<const SpmdVectorSpaceBase<Scalar> >
DefaultSpmdMultiVector<Scalar>::spmdSpace() const
{
#ifdef THYRA_DEFAULT_SPMD_MULTI_VECTOR_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nSpmdMultiVectorStd<Scalar>::spmdSpace() const called!\n";
#endif
  return spmdRangeSpace_;
}


// Overridden protected functions from MultiVectorBase


template<class Scalar>
RCP<VectorBase<Scalar> >
DefaultSpmdMultiVector<Scalar>::nonconstColImpl(Ordinal j)
{
#ifdef THYRA_DEFAULT_SPMD_MULTI_VECTOR_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nSpmdMultiVectorStd<Scalar>::col() called!\n";
#endif
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( !( 0 <= j && j < this->domain()->dim() ) );
#endif
  return Teuchos::rcp(
    new DefaultSpmdVector<Scalar>(
      spmdRangeSpace_,
      localValues_.persistingView(j*leadingDim_,spmdRangeSpace_->localSubDim()),
      1
      )
    );
  //return Teuchos::rcp(new DefaultVectorMultiVector<Scalar>(subView(Range1D(j,j))));
}


template<class Scalar>
RCP<const MultiVectorBase<Scalar> >
DefaultSpmdMultiVector<Scalar>::contigSubViewImpl(
  const Range1D& col_rng_in
  ) const
{
  const Range1D colRng = this->validateColRange(col_rng_in);
  return Teuchos::rcp(
    new DefaultSpmdMultiVector<Scalar>(
      spmdRangeSpace_,
      Teuchos::rcp_dynamic_cast<const ScalarProdVectorSpaceBase<Scalar> >(
        spmdRangeSpace_->smallVecSpcFcty()->createVecSpc(colRng.size())
        ,true
        ),
      localValues_.persistingView(colRng.lbound()*leadingDim_,colRng.size()*spmdRangeSpace_->localSubDim()),
      leadingDim_
      )
    );
}


template<class Scalar>
RCP<MultiVectorBase<Scalar> >
DefaultSpmdMultiVector<Scalar>::nonconstContigSubViewImpl(
  const Range1D& col_rng_in
  )
{
  return Teuchos::rcp_const_cast<MultiVectorBase<Scalar> >(
    this->contigSubViewImpl(col_rng_in));
  // Have the nonconst version call the const version.  Note that in this case
  // we just need to take the const off of the returned MultiVectorBase object
  // because the localValues is already handled as nonconst.  This is the
  // perfect instance where the advice in Item 3 in "Effective C++ 3rd
  // edition" where Scott Meyers recommends having the nonconst version call
  // the const version.
}


template<class Scalar>
RCP<const MultiVectorBase<Scalar> >
DefaultSpmdMultiVector<Scalar>::nonContigSubViewImpl(
  const ArrayView<const int> &cols
  ) const
{
  THYRA_DEBUG_ASSERT_MV_COLS("nonContigSubViewImpl(cols)", cols);
  const int numCols = cols.size();
  const ArrayRCP<Scalar> localValuesView = createContiguousCopy(cols);
  return defaultSpmdMultiVector<Scalar>(
    spmdRangeSpace_,
    createSmallScalarProdVectorSpaceBase<Scalar>(spmdRangeSpace_, numCols),
    localValuesView
    );
}


template<class Scalar>
RCP<MultiVectorBase<Scalar> >
DefaultSpmdMultiVector<Scalar>::nonconstNonContigSubViewImpl(
  const ArrayView<const int> &cols )
{
  THYRA_DEBUG_ASSERT_MV_COLS("nonContigSubViewImpl(cols)", cols);
  const int numCols = cols.size();
  const ArrayRCP<Scalar> localValuesView = createContiguousCopy(cols);
  const Ordinal localSubDim = spmdRangeSpace_->localSubDim();
  RCP<CopyBackSpmdMultiVectorEntries<Scalar> > copyBackView =
    copyBackSpmdMultiVectorEntries<Scalar>(cols, localValuesView.getConst(),
      localSubDim, localValues_.create_weak(), leadingDim_);
  return Teuchos::rcpWithEmbeddedObjPreDestroy(
    new DefaultSpmdMultiVector<Scalar>(
      spmdRangeSpace_,
      createSmallScalarProdVectorSpaceBase<Scalar>(spmdRangeSpace_, numCols),
      localValuesView),
    copyBackView
    );
}


// Overridden protected members from SpmdMultiVectorBase


template<class Scalar>
void DefaultSpmdMultiVector<Scalar>::getNonconstLocalDataImpl(
  const Ptr<ArrayRCP<Scalar> > &localValues, const Ptr<Ordinal> &leadingDim
  )
{
#ifdef THYRA_DEFAULT_SPMD_MULTI_VECTOR_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nSpmdMultiVectorStd<Scalar>::getLocalDataImpl() called!\n";
#endif
  *localValues = localValues_;
  *leadingDim = leadingDim_;
}


template<class Scalar>
void DefaultSpmdMultiVector<Scalar>::getLocalDataImpl(
  const Ptr<ArrayRCP<const Scalar> > &localValues, const Ptr<Ordinal> &leadingDim
  ) const
{
#ifdef THYRA_DEFAULT_SPMD_MULTI_VECTOR_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nSpmdMultiVectorStd<Scalar>::getLocalData() called!\n";
#endif
  *localValues = localValues_;
  *leadingDim = leadingDim_;
}


// private


template<class Scalar>
ArrayRCP<Scalar>
DefaultSpmdMultiVector<Scalar>::createContiguousCopy(
  const ArrayView<const int> &cols ) const
{
  typedef typename ArrayRCP<Scalar>::const_iterator const_itr_t;
  typedef typename ArrayRCP<Scalar>::iterator itr_t;
  const int numCols = cols.size();
  const Ordinal localSubDim = spmdRangeSpace_->localSubDim();
  ArrayRCP<Scalar> localValuesView = Teuchos::arcp<Scalar>(numCols*localSubDim);
  // Copy to contiguous storage column by column
  const const_itr_t lv = localValues_.begin();
  const itr_t lvv = localValuesView.begin();
  for (int k = 0; k < numCols; ++k) {
    const int col_k = cols[k];
    const const_itr_t lv_k = lv + leadingDim_*col_k;
    const itr_t lvv_k = lvv + localSubDim*k;
    std::copy(lv_k, lv_k+localSubDim, lvv_k);
  }
  return localValuesView;
}


} // end namespace Thyra


#endif // THYRA_DEFAULT_SPMD_MULTI_VECTOR_DEF_HPP
