// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
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
//
#ifndef BELOS_MULTI_VEC_TRAITS_HPP
#define BELOS_MULTI_VEC_TRAITS_HPP

#include "BelosMultiVec.hpp"
#ifdef HAVE_TSFCORE
#include "WorkspacePack.hpp"
#include "TSFCoreExplicitMultiVectorView.hpp"
#include "TSFCoreTestingTools.hpp"
#include "TSFCoreMultiVector.hpp"
#include "TSFCoreMultiVectorStdOps.hpp"
#include "TSFCoreVectorSpace.hpp"
#endif

namespace Belos {

template<class TYPE, class MV>
class MultiVecTraits {};

template<class TYPE>
class MultiVecTraits<TYPE,MultiVec<TYPE> >
{
public:
  ///
  static const MultiVec<TYPE>& c(MultiVec<TYPE>& mv) { return mv; } 

  ///
  static RefCountPtr<MultiVec<TYPE> > Clone( const MultiVec<TYPE>& mv, const int numvecs )
  { return rcp( const_cast<MultiVec<TYPE>&>(mv).Clone(numvecs) ); }
  ///
  static RefCountPtr<MultiVec<TYPE> > CloneCopy( const MultiVec<TYPE>& mv )
  { return rcp( const_cast<MultiVec<TYPE>&>(mv).CloneCopy() ); }
  ///
  //static RefCountPtr<MultiVec<TYPE> > CloneCopy( const MultiVec<TYPE>& mv, int index[], int numvecs )
  //{ return rcp( const_cast<MultiVec<TYPE>&>(mv).CloneCopy(index,numvecs) ); }
  ///
  static RefCountPtr<MultiVec<TYPE> > CloneView( MultiVec<TYPE>& mv, int index[], int numvecs )
  { return rcp( mv.CloneView(index,numvecs) ); }
  ///
  static RefCountPtr<const MultiVec<TYPE> > CloneView( const MultiVec<TYPE>& mv, int index[], int numvecs )
  { return rcp( const_cast<MultiVec<TYPE>&>(mv).CloneView(index,numvecs) ); }
  ///
  static int GetVecLength( const MultiVec<TYPE>& mv )
  { return mv.GetVecLength(); }
  ///
  static int GetNumberVecs( const MultiVec<TYPE>& mv )
  { return mv.GetNumberVecs(); }
  ///
  static void MvTimesMatAddMv( TYPE alpha, const MultiVec<TYPE>& A, 
			       const Teuchos::SerialDenseMatrix<int,TYPE>& B, TYPE beta, MultiVec<TYPE>& mv )
    { mv.MvTimesMatAddMv(alpha,const_cast<MultiVec<TYPE>&>(A),const_cast<Teuchos::SerialDenseMatrix<int,TYPE>&>(B),beta); }
  ///
  static void MvAddMv( TYPE alpha, const MultiVec<TYPE>& A, TYPE beta, const MultiVec<TYPE>& B, MultiVec<TYPE>& mv )
    { mv.MvAddMv(alpha,const_cast<MultiVec<TYPE>&>(A),beta,const_cast<MultiVec<TYPE>&>(B)); }
  ///
  static void MvTransMv( const MultiVec<TYPE>& mv, TYPE alpha, const MultiVec<TYPE>& A, Teuchos::SerialDenseMatrix<int,TYPE>& B )
    { const_cast<MultiVec<TYPE>&>(mv).MvTransMv(alpha,const_cast<MultiVec<TYPE>&>(A),B); }
  ///
  static ReturnType MvNorm( const MultiVec<TYPE>& mv, TYPE *normvec, NormType norm_type = TwoNorm )
    { return const_cast<MultiVec<TYPE>&>(mv).MvNorm(normvec,norm_type); }
  ///
  static void SetBlock( const MultiVec<TYPE>& A, int index[], int numvecs, MultiVec<TYPE>& mv )
    { mv.SetBlock(const_cast<MultiVec<TYPE>&>(A),index,numvecs); }
  ///
  static void MvRandom( MultiVec<TYPE>& mv )
    { mv.MvRandom(); }
  ///
  static void MvInit( MultiVec<TYPE>& mv, TYPE alpha = Teuchos::ScalarTraits<TYPE>::zero() )
    { mv.MvInit(alpha); }
  ///
  static void MvPrint( const MultiVec<TYPE>& mv, ostream& os )
    { const_cast<MultiVec<TYPE>&>(mv).MvPrint(os); }
    
};

#ifdef HAVE_TSFCORE

template<class TYPE>
class MultiVecTraits<TYPE,TSFCore::MultiVector<TYPE> >
{
private:
  ///
  typedef Teuchos::ScalarTraits<TYPE> ST;
  ///
  static void copy(
	    const Teuchos::SerialDenseMatrix<int,TYPE>  &sdm
	    ,TSFCore::MultiVector<TYPE>                 *mv
	    )
  {
    TSFCore::ExplicitMutableMultiVectorView<TYPE> mvv(*mv);
    Teuchos::BLAS<int,TYPE> blas;
    if( mvv.leadingDim() == sdm.stride() ) {
      blas.COPY( sdm.numRows() * sdm.numCols(), &sdm(0,0), 1, &mvv(1,1), 1 );
    }
    else {
      for( int j = 0; j < sdm.numCols(); ++j ) {
	blas.COPY( sdm.numRows(), &sdm(0,j), 1, &mvv(1,j+1), 1 );
      }
    }
  }
  ///
  static void copy(
	    const TSFCore::MultiVector<TYPE>       &mv
	    ,Teuchos::SerialDenseMatrix<int,TYPE>  *sdm_out
	    )
  {
    Teuchos::SerialDenseMatrix<int,TYPE> &sdm = *sdm_out;
    TSFCore::ExplicitMultiVectorView<TYPE> mvv(mv);
    Teuchos::BLAS<int,TYPE> blas;
    if( mvv.leadingDim() == sdm.stride() ) {
      blas.COPY( sdm.numRows() * sdm.numCols(), &mvv(1,1), 1, &sdm(0,0), 1 );
    }
    else {
      for( int j = 0; j < sdm.numCols(); ++j ) {
	blas.COPY( sdm.numRows(), &mvv(1,j+1), 1, &sdm(0,j), 1 );
      }
    }
  }
  ///
  static RefCountPtr<TSFCore::MultiVector<TYPE> > convert(
	  const TSFCore::VectorSpace<TYPE>             &space
	  ,const Teuchos::SerialDenseMatrix<int,TYPE>  &sdm
	  )
  {
    TEST_FOR_EXCEPT( space.dim() != sdm.numRows() );
    RefCountPtr<TSFCore::MultiVector<TYPE> > mv = space.createMembers(sdm.numCols());
    copy( sdm, &*mv );
    return mv;
  }
public:
  ///
  static const TSFCore::MultiVector<TYPE>& c(TSFCore::MultiVector<TYPE>& mv) { return mv; } 
  ///
  static RefCountPtr<TSFCore::MultiVector<TYPE> > Clone( const TSFCore::MultiVector<TYPE>& mv, const int numvecs )
  { return mv.range()->createMembers(numvecs); }
  ///
  static RefCountPtr<TSFCore::MultiVector<TYPE> > CloneCopy( const TSFCore::MultiVector<TYPE>& mv )
  { RefCountPtr<TSFCore::MultiVector<TYPE> > mv_out = Clone(mv,GetNumberVecs(mv)); assign(&*mv_out,mv); return mv_out; }
  ///
  static RefCountPtr<TSFCore::MultiVector<TYPE> > CloneView( TSFCore::MultiVector<TYPE>& mv, int index[], int numvecs )
  {
    namespace wsp = WorkspacePack;
    wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();
    wsp::Workspace<int> index_1(wss,numvecs);
    for( int k = 0; k < numvecs; ++k ) index_1[k] = index[k]+1;
    return mv.subView(numvecs,&index_1[0]);
  }
  ///
  static RefCountPtr<const TSFCore::MultiVector<TYPE> > CloneView( const TSFCore::MultiVector<TYPE>& mv, int index[], int numvecs )
  {
    namespace wsp = WorkspacePack;
    wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();
    wsp::Workspace<int> index_1(wss,numvecs);
    for( int k = 0; k < numvecs; ++k ) index_1[k] = index[k]+1;
    return mv.subView(numvecs,&index_1[0]);
  }
  ///
  static int GetVecLength( const TSFCore::MultiVector<TYPE>& mv )
  { return mv.range()->dim(); }
  ///
  static int GetNumberVecs( const TSFCore::MultiVector<TYPE>& mv )
  { return mv.domain()->dim(); }
  ///
  static void MvTimesMatAddMv( TYPE alpha, const TSFCore::MultiVector<TYPE>& A, 
			       const Teuchos::SerialDenseMatrix<int,TYPE>& B, TYPE beta, TSFCore::MultiVector<TYPE>& mv )
    {
      // mv = beta * mv + alpha * A * B
      A.apply( TSFCore::NOTRANS, *convert(*A.domain(),B), &mv, alpha, beta );
    }
  ///
  static void MvAddMv( TYPE alpha, const TSFCore::MultiVector<TYPE>& A, TYPE beta, const TSFCore::MultiVector<TYPE>& B, TSFCore::MultiVector<TYPE>& mv )
    {
      // mv = alpha*A + beta*B
      if ( &A == &mv ) {
	//
	// *this *= alpha
	TSFCore::scale( alpha, &mv );
	//
	// *this += beta * B
	TSFCore::update( beta, B, &mv );
	//
      } else if ( &B == &mv ) { 
	//
	// *this *=beta
	TSFCore::scale( beta, &mv );
	//
	// *this += alpha * A
	TSFCore::update( alpha, A, &mv );
	//
      } else {
	// *this <- A
	TSFCore::assign( &mv, A );
	//
	// *this *= alpha
	TSFCore::scale( alpha, &mv );
	//
	// *this += beta * B
	TSFCore::update( beta, B, &mv );
      }
    }
  ///
  static void MvTransMv( const TSFCore::MultiVector<TYPE>& mv, TYPE alpha, const TSFCore::MultiVector<TYPE>& A, Teuchos::SerialDenseMatrix<int,TYPE>& B )
    {
      // B = alpha * A' * mv
      RefCountPtr<TSFCore::MultiVector<TYPE> > B_mv = mv.domain()->createMembers(A.domain()->dim());
      A.apply( TSFCore::TRANS, mv, &*B_mv, alpha, ST::zero() );
      copy( *B_mv, &B );
    }
  ///
  static ReturnType MvNorm( const TSFCore::MultiVector<TYPE>& mv, TYPE *normvec, NormType norm_type = TwoNorm )
  {
    const int numVecs = GetNumberVecs(mv);
    if (normvec) {
      switch( norm_type ) {
      case ( OneNorm ) :
	for (int i=0; i<numVecs; i++)
	  normvec[i] = TSFCore::norm_1( *(mv.col(i+1)) );
	return Ok;
      case ( TwoNorm ) :
	for (int i=0; i<numVecs; i++)
	  normvec[i] = TSFCore::norm_2( *(mv.col(i+1)) );
	return Ok;
      case ( InfNorm ) :
	for (int i=0; i<numVecs; i++) 
	  normvec[i] = TSFCore::norm_inf( *(mv.col(i+1)) );
	return Ok;
      default :
	return Undefined;
      }
    }
    return Undefined;
  }
  ///
  static void SetBlock( const TSFCore::MultiVector<TYPE>& A, int index[], int numvecs, TSFCore::MultiVector<TYPE>& mv )
  {
    for( int k = 0; k < numvecs; ++k )
      assign( &*mv.col(k+1), *A.col(index[k]+1) );
  }
  ///
  static void MvRandom( TSFCore::MultiVector<TYPE>& mv )
  {
    TSFCore::randomize( 0.0, 1.0, &mv );
  }
  ///
  static void MvInit( TSFCore::MultiVector<TYPE>& mv, TYPE alpha = Teuchos::ScalarTraits<TYPE>::zero() )
  {
    TSFCore::assign( &mv, alpha );
  }
  ///
  static void MvPrint( const TSFCore::MultiVector<TYPE>& mv, ostream& os )
  {
    os << mv;
  }
};

#endif // HAVE_TSFCORE

} // namespace Belos

#endif // BELOS_MULTI_VEC_TRAITS_HPP
