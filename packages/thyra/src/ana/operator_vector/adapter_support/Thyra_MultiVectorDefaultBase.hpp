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

#ifndef THYRA_MULTI_VECTOR_DEFAULT_BASE_HPP
#define THYRA_MULTI_VECTOR_DEFAULT_BASE_HPP

#include "Thyra_MultiVectorDefaultBaseDecl.hpp"
#include "Thyra_LinearOpDefaultBase.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_AssertOp.hpp"
#include "Thyra_MultiVectorCols.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_TestForException.hpp"

namespace Thyra {


// Sub-view methods

template<class Scalar>
Teuchos::RefCountPtr<const MultiVectorBase<Scalar> >
MultiVectorDefaultBase<Scalar>::subView( const Range1D& colRng_in ) const
{
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore        *wss      = Teuchos::get_default_workspace_store().get();
  const VectorSpaceBase<Scalar>  &domain   = *this->domain();
  const VectorSpaceBase<Scalar>  &range    = *this->range();
  const Index                    dimDomain = domain.dim();
  const Range1D                  colRng    = RangePack::full_range(colRng_in,1,dimDomain);
  if( colRng.lbound() == 1 && static_cast<Index>(colRng.ubound()) == dimDomain )
    return Teuchos::rcp(this,false); // Takes all of the columns!
  if( colRng.size() ) {
    // We have to create a view of a subset of the columns
    Workspace< Teuchos::RefCountPtr< VectorBase<Scalar> > >  col_vecs(wss,colRng.size());
    for( Index j = colRng.lbound(); j <= colRng.ubound(); ++j )
      col_vecs[j-colRng.lbound()] = Teuchos::rcp_const_cast<VectorBase<Scalar> >(this->col(j));
    return Teuchos::rcp(new MultiVectorCols<Scalar>(this->range(),range.smallVecSpcFcty()->createVecSpc(colRng.size()),&col_vecs[0]));
  }
  return Teuchos::null; // There was an empty set in colRng_in!
}

template<class Scalar>
Teuchos::RefCountPtr<MultiVectorBase<Scalar> >
MultiVectorDefaultBase<Scalar>::subView( const Range1D& colRng_in )
{
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore        *wss      = Teuchos::get_default_workspace_store().get();
  const VectorSpaceBase<Scalar>  &domain   = *this->domain();
  const VectorSpaceBase<Scalar>  &range    = *this->range();
  const Index                    dimDomain = domain.dim();
  const Range1D                  colRng    = RangePack::full_range(colRng_in,1,dimDomain);
  if( colRng.lbound() == 1 && static_cast<Index>(colRng.ubound()) == dimDomain )
    return Teuchos::rcp(this,false); // Takes all of the columns!
  if( colRng.size() ) {
    // We have to create a view of a subset of the columns
    Workspace< Teuchos::RefCountPtr< VectorBase<Scalar> > >  col_vecs(wss,colRng.size());
    for( Index j = colRng.lbound(); j <= colRng.ubound(); ++j )
      col_vecs[j-colRng.lbound()] = this->col(j);
    return Teuchos::rcp(new MultiVectorCols<Scalar>(this->range(),range.smallVecSpcFcty()->createVecSpc(colRng.size()),&col_vecs[0]));
  }
  return Teuchos::null; // There was an empty set in colRng_in!
}

template<class Scalar>
Teuchos::RefCountPtr<const MultiVectorBase<Scalar> >
MultiVectorDefaultBase<Scalar>::subView( const int numCols, const int cols[] ) const
{
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore        *wss      = Teuchos::get_default_workspace_store().get();
  const VectorSpaceBase<Scalar>  &range    = *this->range();
#ifdef _DEBUG
  const VectorSpaceBase<Scalar>  &domain   = *this->domain();
  const Index                    dimDomain = domain.dim();
  const char msg_err[] = "MultiVectorDefaultBase<Scalar>::subView(numCols,cols[]): Error!";
   TEST_FOR_EXCEPTION( numCols < 1 || dimDomain < numCols, std::invalid_argument, msg_err );
#endif
  // We have to create a view of a subset of the columns
  Workspace< Teuchos::RefCountPtr< VectorBase<Scalar> > > col_vecs(wss,numCols);
  for( int k = 0; k < numCols; ++k ) {
    const int col_k = cols[k];
#ifdef _DEBUG
    TEST_FOR_EXCEPTION(
      col_k < 1 || dimDomain < col_k, std::invalid_argument
      ,msg_err << " col["<<k<<"] = " << col_k << " is not in the range [1,"<<dimDomain<<"]!"
      );
#endif
    col_vecs[k] = Teuchos::rcp_const_cast<VectorBase<Scalar> >(this->col(col_k));
  }
  return Teuchos::rcp(new MultiVectorCols<Scalar>(this->range(),range.smallVecSpcFcty()->createVecSpc(numCols),&col_vecs[0]));
}

template<class Scalar>
Teuchos::RefCountPtr<MultiVectorBase<Scalar> >
MultiVectorDefaultBase<Scalar>::subView( const int numCols, const int cols[] )
{
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore        *wss      = Teuchos::get_default_workspace_store().get();
  const VectorSpaceBase<Scalar>  &range    = *this->range();
#ifdef _DEBUG
  const VectorSpaceBase<Scalar>  &domain   = *this->domain();
  const Index                    dimDomain = domain.dim();
  const char msg_err[] = "MultiVectorDefaultBase<Scalar>::subView(numCols,cols[]): Error!";
   TEST_FOR_EXCEPTION( numCols < 1 || dimDomain < numCols, std::invalid_argument, msg_err );
#endif
  // We have to create a view of a subset of the columns
  Workspace< Teuchos::RefCountPtr< VectorBase<Scalar> > > col_vecs(wss,numCols);
  for( int k = 0; k < numCols; ++k ) {
    const int col_k = cols[k];
#ifdef _DEBUG
    TEST_FOR_EXCEPTION(
      col_k < 1 || dimDomain < col_k, std::invalid_argument
      ,msg_err << " col["<<k<<"] = " << col_k << " is not in the range [1,"<<dimDomain<<"]!"
      );
#endif
    col_vecs[k] = this->col(col_k);
  }
  return Teuchos::rcp(new MultiVectorCols<Scalar>(this->range(),range.smallVecSpcFcty()->createVecSpc(numCols),&col_vecs[0]));
}

} // end namespace Thyra

#endif // THYRA_MULTI_VECTOR_DEFAULT_BASE_HPP
