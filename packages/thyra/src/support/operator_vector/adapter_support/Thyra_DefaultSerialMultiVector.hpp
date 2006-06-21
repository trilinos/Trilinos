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

#ifndef THYRA_SERIAL_MULTI_VECTOR_STD_HPP
#define THYRA_SERIAL_MULTI_VECTOR_STD_HPP

// Define to make some verbose output
//#define THYRA_SERIAL_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT

#include "Thyra_DefaultSerialMultiVectorDecl.hpp"
#include "Thyra_SerialMultiVectorBase.hpp"
#include "Thyra_DefaultSerialVector.hpp"

namespace Thyra {

//
// Simple utility class to copy back multi-vector entries
//

template<class Scalar>
class CopyBackSerialMultiVectorEntries {
public:
  CopyBackSerialMultiVectorEntries(
    const int numCols, const int cols[]
    ,const Scalar *valuesView, const Index dim
    ,Scalar *values,const Index leadingDim
    )
    :cols_(cols,cols+numCols)
    ,valuesView_(valuesView),dim_(dim)
    ,values_(values),leadingDim_(leadingDim)
    {}
  ~CopyBackSerialMultiVectorEntries()
    {
      //std::cout << "\nEntering CopyBackSerialMultiVectorEntries::~CopyBackSerialMultiVectorEntries()...\n";
      // Copy from contiguous storage
      const int    numCols = cols_.size();
      const Scalar *lvv    = &*valuesView_;
      Scalar       *lv     = &*values_;
      for( int k = 0; k < numCols; ++k ) {
        const int col_k = cols_[k];
        const Scalar *lvv_k  = lvv + dim_*k;
        Scalar       *lv_k   = lv + leadingDim_*col_k;
        //std::cout << "\nlvv_k = ["; for( int j = 0; j < dim_; ++j ) std::cout << lvv_k[j] << ","; std::cout << "]\n";
        std::copy( lvv_k, lvv_k + dim_, lv_k );
        //std::cout << "\nlv_k = ["; for( int j = 0; j < dim_; ++j ) std::cout << lv_k[j] << ","; std::cout << "]\n";
      }
      //std::cout << "\nLeaving CopyBackSerialMultiVectorEntries::~CopyBackSerialMultiVectorEntries()\n";
    }
private:
  std::vector<int>   cols_;
  const Scalar       *valuesView_;
  Index              dim_;
  Scalar             *values_;
  Index              leadingDim_;
  // Not defined and not to be called
  CopyBackSerialMultiVectorEntries();
  CopyBackSerialMultiVectorEntries(const CopyBackSerialMultiVectorEntries&);
  CopyBackSerialMultiVectorEntries& operator=(const CopyBackSerialMultiVectorEntries&);
};

//
// DefaultSerialMultiVector
//

// Constructors/initializers/accessors

template<class Scalar>
DefaultSerialMultiVector<Scalar>::DefaultSerialMultiVector()
  :leadingDim_(0)
{}

template<class Scalar>
DefaultSerialMultiVector<Scalar>::DefaultSerialMultiVector(
  const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >   &range
  ,const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >  &domain
  )
{
  initialize(range,domain);
}

template<class Scalar>
DefaultSerialMultiVector<Scalar>::DefaultSerialMultiVector(
  const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >   &range
  ,const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >  &domain
  ,const Teuchos::RefCountPtr<Scalar>                                    &values
  ,const Index                                                           leadingDim
  )
{
  initialize(range,domain,values,leadingDim);
}

template<class Scalar>
void DefaultSerialMultiVector<Scalar>::initialize(
  const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >   &range
  ,const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >  &domain
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(range.get()==NULL);
  TEST_FOR_EXCEPT(domain.get()==NULL);
#endif
  const Index numRows = range->dim(), numCols = domain->dim();
  initialize(
    range, domain
    ,Teuchos::rcp(new Scalar[numRows*numCols],Teuchos::DeallocArrayDelete<Scalar>(),true)
    ,numRows
    );
}

template<class Scalar>
void DefaultSerialMultiVector<Scalar>::initialize(
  const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >   &range
  ,const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >  &domain
  ,const Teuchos::RefCountPtr<Scalar>                                    &values
  ,const Index                                                           leadingDim
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(range.get()==NULL);
  TEST_FOR_EXCEPT(domain.get()==NULL);
  TEST_FOR_EXCEPT(values.get()==NULL);
  TEST_FOR_EXCEPT(leadingDim < range->dim());
#endif
  range_       = range;
  domain_      = domain;
  values_      = values;
  leadingDim_  = leadingDim;
  numRows_     = range->dim();
  numCols_     = domain->dim();
  this->updateSpace();
}

template<class Scalar>
void DefaultSerialMultiVector<Scalar>::uninitialize(
  Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >        *range
  ,Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >       *domain
  ,Teuchos::RefCountPtr<Scalar>                                         *values
  ,Index                                                                *leadingDim
  )
{
  if(range)         *range         = range_;
  if(domain)        *domain        = domain_;
  if(values)        *values        = values_;
  if(leadingDim)    *leadingDim    = leadingDim_;

  range_       = Teuchos::null;
  domain_      = Teuchos::null;
  values_      = Teuchos::null;
  leadingDim_  = 0;

  this->updateSpace();
}

// Overridden from Teuchos::Describable

template<class Scalar>
std::string DefaultSerialMultiVector<Scalar>::description() const
{
  return (std::string("DefaultSerialMultiVector<") + Teuchos::ScalarTraits<Scalar>::name() + std::string(">"));
}

// Overridden from EuclideanLinearOpBase

template<class Scalar>
Teuchos::RefCountPtr< const ScalarProdVectorSpaceBase<Scalar> >
DefaultSerialMultiVector<Scalar>::rangeScalarProdVecSpc() const
{
#ifdef THYRA_SERIAL_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nSerialMultiVectorStd<Scalar>::range() const called!\n";
#endif
  return range_;
}

template<class Scalar>
Teuchos::RefCountPtr< const ScalarProdVectorSpaceBase<Scalar> >
DefaultSerialMultiVector<Scalar>::domainScalarProdVecSpc() const
{
#ifdef THYRA_SERIAL_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nSerialMultiVectorStd<Scalar>::domain() const called!\n";
#endif
  return domain_;
}

// Overridden from MultiVectorBase

template<class Scalar>
Teuchos::RefCountPtr<VectorBase<Scalar> >
DefaultSerialMultiVector<Scalar>::col(Index j)
{
  using Teuchos::rcp;
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( !( 0 <= j && j < numCols_ ) );
#endif
#ifdef THYRA_SERIAL_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nSerialMultiVectorStd<Scalar>::col() called!\n";
#endif
  return rcp(new DefaultSerialVector<Scalar>(rcp((&*values_)+j*leadingDim_,false),1,numRows_,range_));
}

template<class Scalar>
Teuchos::RefCountPtr<MultiVectorBase<Scalar> >
DefaultSerialMultiVector<Scalar>::subView( const Range1D& col_rng_in )
{
#ifdef THYRA_SERIAL_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nSerialMultiVectorStd<Scalar>::subView() called!\n";
#endif
  const Range1D colRng = this->validateColRange(col_rng_in);
  return Teuchos::rcp(
    new DefaultSerialMultiVector<Scalar>(
      range_
      ,Teuchos::rcp_dynamic_cast<const ScalarProdVectorSpaceBase<Scalar> >(range_->smallVecSpcFcty()->createVecSpc(colRng.size()),true)
      ,Teuchos::rcp( (&*values_) + colRng.lbound()*leadingDim_, false )
      ,leadingDim_
      )
    );
}

template<class Scalar>
Teuchos::RefCountPtr<const MultiVectorBase<Scalar> >
DefaultSerialMultiVector<Scalar>::subView( const int numCols, const int cols[] ) const
{
  Teuchos::RefCountPtr<Scalar> valuesView = createContiguousCopy(numCols,cols);
  return 
    Teuchos::rcp(
      new DefaultSerialMultiVector<Scalar>(
        range_
        ,Teuchos::rcp_dynamic_cast<const ScalarProdVectorSpaceBase<Scalar> >(
          range_->smallVecSpcFcty()->createVecSpc(numCols)
          ,true)
        ,valuesView
        ,range_->dim()
        )
      );
}

template<class Scalar>
Teuchos::RefCountPtr<MultiVectorBase<Scalar> >
DefaultSerialMultiVector<Scalar>::subView( const int numCols, const int cols[] )
{
  Teuchos::RefCountPtr<Scalar> valuesView = createContiguousCopy(numCols,cols);
  const Index dim = range_->dim();
  Teuchos::RefCountPtr<MultiVectorBase<Scalar> >
    view = Teuchos::rcp(
      new DefaultSerialMultiVector<Scalar>(
        range_
        ,Teuchos::rcp_dynamic_cast<const ScalarProdVectorSpaceBase<Scalar> >(
          range_->smallVecSpcFcty()->createVecSpc(numCols)
          ,true)
        ,valuesView
        ,dim
        )
      );
  Teuchos::RefCountPtr<CopyBackSerialMultiVectorEntries<Scalar> >
    copyBackView
    = Teuchos::rcp(
      new CopyBackSerialMultiVectorEntries<Scalar>(
        numCols,cols,&*valuesView,dim,&*values_,leadingDim_
        )
      );
  Teuchos::set_extra_data(copyBackView,"copyBackView",&view,Teuchos::PRE_DESTROY);
  return view;
}

// Overridden from SerialMultiVectorBase

template<class Scalar>
void DefaultSerialMultiVector<Scalar>::getData( const Scalar **values, Index *leadingDim ) const
{
#ifdef THYRA_SERIAL_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nSerialMultiVectorStd<Scalar>::getData() const called!\n";
#endif
  *values     = &*values_;
  *leadingDim = leadingDim_;
}

template<class Scalar>
void DefaultSerialMultiVector<Scalar>::freeData( const Scalar *values ) const
{
#ifdef THYRA_SERIAL_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nSerialMultiVectorStd<Scalar>::freeData() called!\n";
#endif
}

template<class Scalar>
void DefaultSerialMultiVector<Scalar>::getData( Scalar **values, Index *leadingDim )
{
#ifdef THYRA_SERIAL_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nSerialMultiVectorStd<Scalar>::getData() called!\n";
#endif
  *values     = &*values_;
  *leadingDim = leadingDim_;
}

template<class Scalar>
void DefaultSerialMultiVector<Scalar>::commitData( Scalar *values )
{
#ifdef THYRA_SERIAL_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nSerialMultiVectorStd<Scalar>::commitData() called!\n";
#endif
}

// private

template<class Scalar>
Teuchos::RefCountPtr<Scalar>
DefaultSerialMultiVector<Scalar>::createContiguousCopy( const int numCols, const int cols[] ) const
{
#ifdef TEUCHOS_DEBUG
  const VectorSpaceBase<Scalar>  &domain   = *domain_;
  const Index                    dimDomain = domain.dim();
  const char msg_err[] = "MultiVectorDefaultBase<Scalar>::subView(numCols,cols[]): Error!";
  TEST_FOR_EXCEPTION( numCols < 1 || dimDomain < numCols, std::invalid_argument, msg_err );
#endif
  const Index dim = range_->dim();
  Teuchos::RefCountPtr<Scalar>
    valuesView = Teuchos::rcp( new Scalar[numCols*dim], Teuchos::DeallocArrayDelete<Scalar>(), true );
  // Copy to contiguous storage
  const Scalar *lv   = &*values_;
  Scalar       *lvv  = &*valuesView;
  for( int k = 0; k < numCols; ++k ) {
    const int col_k = cols[k];
#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPTION(
      col_k < 1 || dimDomain < col_k, std::invalid_argument
      ,msg_err << " col["<<k<<"] = " << col_k << " is not in the range [1,"<<dimDomain<<"]!"
      );
#endif
    const Scalar *lv_k   = lv + leadingDim_*col_k;
    Scalar       *lvv_k  = lvv + dim*k;
    std::copy( lv_k, lv_k + dim, lvv_k );
  }
  return valuesView;
}

} // end namespace Thyra

#endif // THYRA_SERIAL_MULTI_VECTOR_STD_HPP
