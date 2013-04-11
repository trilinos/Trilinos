/*
 * KokkosArray_CRSMatrix.h
 *
 *  Created on: Jul 30, 2012
 *      Author: crtrott
 */

#ifndef KOKKOSARRAY_CRSMATRIX_H_
#define KOKKOSARRAY_CRSMATRIX_H_

#include <KokkosArray_View.hpp>
#include <KokkosArray_Cuda.hpp>
#include <KokkosArray_Macros.hpp>
#include <KokkosArray_CrsArray.hpp>
#include <KokkosArray_MultiVector.hpp>
#ifdef KOKKOS_USE_CUSPARSE
#include <cusparse_v2.h>
#include <Kokkos_CRSMatrix_CuSparse.hpp>
#endif
namespace KokkosArray {

//ToDo: Check Type compatibility for Kernel Calls

#define UNROLL_LIMIT 16
template< typename ScalarType , typename OrdinalType, class Device >
class CrsMatrix {
public:
  typedef Device      device_type ;
  typedef ScalarType  scalar_type ;
  typedef OrdinalType  ordinal_type ;

  typedef CrsMatrix<ScalarType, OrdinalType,KokkosArray::Host> HostMirror;
  typedef KokkosArray::CrsArray<OrdinalType, KokkosArray::LayoutLeft, Device> CrsArrayType;
  typedef typename CrsArrayType::entries_type   index_type;
  typedef typename CrsArrayType::row_map_type   row_map_type;
  typedef KokkosArray::View< scalar_type* , device_type >   values_type ;

  ordinal_type numRows;
  ordinal_type numCols;
  ordinal_type nnz;
  ordinal_type nentries;

#ifdef KOKKOS_USE_CUSPARSE
  cusparseHandle_t cusparse_handle;
  cusparseMatDescr_t cusparse_descr;
#endif
  CrsArrayType graph;
  values_type  values;

  CrsMatrix() {numRows = 0; numCols = 0; nnz = 0; allocated=0;};
  CrsMatrix(const std::string &label, OrdinalType nrows, OrdinalType ncols, OrdinalType annz)
    {
	  allocated = 0;
	  create(label, nrows,ncols,annz);
#ifdef KOKKOS_USE_CUSPARSE
cusparseCreate(&cusparse_handle);
cusparseCreateMatDescr(&cusparse_descr);
#endif
    };
  CrsMatrix(const std::string &label, OrdinalType nrows, OrdinalType ncols, OrdinalType annz, ScalarType* val, OrdinalType* rows, OrdinalType* cols, bool pad=false)
    {
	  allocated = 0;
	  if(pad)
        import_padded(label,nrows,ncols,annz,val,rows,cols);
	  else
	    import(label,nrows,ncols,annz,val,rows,cols);

      #ifdef KOKKOS_USE_CUSPARSE
	  cusparseCreate(&cusparse_handle);
	  cusparseCreateMatDescr(&cusparse_descr);
      #endif
    };

  void create(const std::string &label, OrdinalType nrows, OrdinalType ncols, OrdinalType annz);
  void import(const std::string &label, OrdinalType nrows, OrdinalType ncols, OrdinalType annz, ScalarType* val, OrdinalType* rows, OrdinalType* cols);
  void import_padded(const std::string &label, OrdinalType nrows, OrdinalType ncols, OrdinalType annz, ScalarType* val, OrdinalType* rows, OrdinalType* cols);
  void generate(const std::string &label, OrdinalType nrows, OrdinalType ncols,OrdinalType target_nnz,OrdinalType varianz_nel_row, OrdinalType width_row);
private:
  OrdinalType allocated;


};

template< typename ScalarType , typename OrdinalType, class Device >
void CrsMatrix<ScalarType ,OrdinalType, Device >::create(const std::string &label, OrdinalType nrows, OrdinalType ncols, OrdinalType annz)
{
  values = values_type("CrsMatrix::values",annz);
  //graph.entries = index_type("CrsMatrix::Graph",annz);
  //graph.row_map = row_map_type("CrsMatrix::rowPtr",nrows+1);
  numRows = nrows;
  numCols = ncols;
  nnz = annz;
  nentries = annz;
  allocated = 1;
}

template< typename ScalarType , typename OrdinalType, class Device >
void CrsMatrix<ScalarType ,OrdinalType, Device >::import(const std::string &label, OrdinalType nrows, OrdinalType ncols, OrdinalType annz, ScalarType* val, OrdinalType* rows, OrdinalType* cols)
{
  std::string str = label;
  values = values_type(str.append(".values"),annz);

  numRows = nrows;
  numCols = ncols;
  nnz = annz;
  nentries = annz;

  std::vector<int> row_lengths (numRows,0);
  for(int i = 0; i<numRows ; i++)
	  row_lengths[i] = rows[i+1] - rows[i];

  str = label;
  graph = KokkosArray::create_crsarray<CrsArrayType>(str.append(".graph"),row_lengths);
  typename values_type::HostMirror h_values = KokkosArray::create_mirror_view(values);
  typename index_type::HostMirror h_entries = KokkosArray::create_mirror_view(graph.entries);

  for(OrdinalType i=0; i<nnz;i++)
  {
	h_values(i) = val[i];
	h_entries(i) = cols[i];
  }
  KokkosArray::deep_copy(values,h_values);
  KokkosArray::deep_copy(graph.entries,h_entries);
  allocated =1;
}

template< typename ScalarType , typename OrdinalType, class Device >
void CrsMatrix<ScalarType ,OrdinalType, Device >::import_padded(const std::string &label,OrdinalType nrows, OrdinalType ncols, OrdinalType annz, ScalarType* val, OrdinalType* rows, OrdinalType* cols)
{

 /* size_t align = 8;//Device::memory_space::preferred_alignment(sizeof(ScalarType),1);

  rowPtr = ordinals_type("CrsMatrix::rowPtr",nrows);
  rowEndPtr = ordinals_type("CrsMatrix::rowEndPtr",nrows);
  typename ordinals_type::HostMirror h_rowPtr = KokkosArray::create_mirror_view(rowPtr);
  typename ordinals_type::HostMirror h_rowEndPtr = KokkosArray::create_mirror_view(rowEndPtr);
  size_t count=0;
  for(OrdinalType i=0; i<nrows;i++)
  {
	size_t elementcount = rows[i+1]-rows[i];
	h_rowPtr(i) = count;
	h_rowEndPtr(i) = count + elementcount;
	count += ((elementcount+align-1)/align)*align;
  }

  nentries = count;
  nnz = annz;
  numRows = nrows;
  numCols = ncols;
  values = values_type("CrsMatrix::values",nentries);
  graph.entries = ordinals_type("CrsMatrix::graph.entries",nentries);
  typename values_type::HostMirror h_values = KokkosArray::create_mirror_view(values);
  typename ordinals_type::HostMirror h_graph.entries = KokkosArray::create_mirror_view(graph.entries);

  for(OrdinalType row=0; row<numRows;row++)
  {
	int j=h_rowPtr(row);
	int k=rows[row];
	for(; j<h_rowEndPtr(row);j++,k++)
	{
	  h_values(j) = val[k];
	  h_colInd(j) = cols[k];
	}
  }
  KokkosArray::deep_copy(values,h_values);
  KokkosArray::deep_copy(colInd,h_colInd);
  KokkosArray::deep_copy(rowPtr,h_rowPtr);
  KokkosArray::deep_copy(rowEndPtr,h_rowEndPtr);*/

}

template< typename ScalarType , typename OrdinalType, class Device >
void CrsMatrix<ScalarType ,OrdinalType, Device >::generate(const std::string &label, OrdinalType nrows, OrdinalType ncols, OrdinalType target_nnz, OrdinalType varianz_nel_row, OrdinalType width_row)
{
  graph.row_map = row_map_type("CrsMatrix::rowPtr",nrows+1);
  typename row_map_type::HostMirror h_row_map = KokkosArray::create_mirror_view(graph.row_map);

  OrdinalType elements_per_row = target_nnz/nrows;
  srand(13721);
  h_row_map(0) = 0;
  for(int row=0;row<nrows;row++)
  {
    int varianz = (1.0*rand()/INT_MAX-0.5)*varianz_nel_row;
    h_row_map(row+1) = h_row_map(row) + elements_per_row+varianz;
  }
  nnz = h_row_map(nrows);
  values = values_type("CrsMatrix::values",nnz);
  graph.entries = index_type("CrsMatrix::colInd",nnz);
  typename values_type::HostMirror h_values = KokkosArray::create_mirror_view(values);
  typename index_type::HostMirror h_entries = KokkosArray::create_mirror_view(graph.entries);
  for(int row=0;row<nrows;row++)
  {
	 for(int k=h_row_map(row);k<h_row_map(row+1);row++)
	 {
		int pos = (1.0*rand()/INT_MAX-0.5)*width_row;
		if(pos<0) pos+=ncols;
		if(pos>=ncols) pos-=ncols;
		h_entries(k) = pos;
		h_values(k) = 100.0*rand()/INT_MAX-50.0;
	 }
  }
  KokkosArray::deep_copy(values,h_values);
  KokkosArray::deep_copy(graph.entries,h_entries);
  KokkosArray::deep_copy(graph.row_map,h_row_map);
}

template<class RangeVector,class CrsMatrix,class DomainVector,class CoeffVector1,class CoeffVector2, int doalpha = 1,int dobeta = 0>
struct MV_MultiplyFunctor
{
  typedef typename CrsMatrix::device_type                   device_type ;
  typedef typename CrsMatrix::ordinal_type                    size_type ;
  typedef typename CrsMatrix::scalar_type                    scalar_type ;
  typedef typename KokkosArray::View<scalar_type* ,typename CrsMatrix::device_type > range_values;

  CoeffVector1 beta;
  CoeffVector2 alpha;
  CrsMatrix  m_A ;
  DomainVector  m_x ;
  RangeVector  m_y ;
  size_type n;

  //--------------------------------------------------------------------------

  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const size_type iRow ) const
  {
	scalar_type sum[16];


	//do multiples of 16
	size_type kk = 0;
	for ( ; kk < n-16 ; kk+=16 ) {
    #pragma vector aligned
    #pragma ivdep
    for ( size_type k = 0 ; k < 16 ; ++k )
	  sum[k] =  0;

    if(doalpha!=-1)
    {
      const size_type iEntryBegin = m_A.graph.row_map(iRow);
      const size_type iEntryEnd   = m_A.graph.row_map(iRow+1);
      for ( size_type iEntry = iEntryBegin ; iEntry < iEntryEnd ; ++iEntry ) {
    	const scalar_type val = m_A.values(iEntry);
    	const size_type ind = m_A.graph.entries(iEntry);
        #pragma vector aligned
        #pragma ivdep
        #pragma unroll
        for ( size_type k = 0 ; k < 16 ; k++ )
    	   sum[k] +=  val*m_x(ind,kk+k);
      }
    }
    else
    {
      const size_type iEntryBegin = m_A.graph.row_map(iRow);
      const size_type iEntryEnd   = m_A.graph.row_map(iRow+1);
      for ( size_type iEntry = iEntryBegin ; iEntry < iEntryEnd ; ++iEntry ) {
    	const scalar_type val = m_A.values(iEntry);
    	const size_type ind = m_A.graph.entries(iEntry);
        #pragma vector aligned
        #pragma ivdep
        for ( size_type k = 0 ; k < 16 ; ++k )
    	  sum[k] -=  val*m_x(ind,kk+k);
      }
    }

    if(doalpha*doalpha!=1) {
      #pragma vector aligned
      #pragma ivdep
      for ( size_type k = 0 ; k < 16 ; ++k )
        sum[k] *= alpha(kk+k);
    }

    if(dobeta == 0) {
      #pragma vector aligned
      #pragma ivdep
      for ( size_type k = 0 ; k < 16 ; ++k )
        m_y(iRow,kk+k) = sum[k] ;
    } else if(dobeta == 1) {
      #pragma vector aligned
      #pragma ivdep
      for ( size_type k = 0 ; k < 16 ; ++k )
        m_y(iRow,kk+k) += sum[k] ;
    } else if(dobeta == -1) {
      #pragma vector aligned
      #pragma ivdep
      for ( size_type k = 0 ; k < 16 ; ++k )
        m_y(iRow,kk+k) = -m_y(iRow,kk+k) +  sum[k] ;
    } else {
      #pragma vector aligned
      #pragma ivdep
      for ( size_type k = 0 ; k < 16 ; ++k )
        m_y(iRow,kk+k) = beta(kk+k)*m_y(iRow,kk+k) + sum[k] ;
    }
  }

    //do remainder
	const int remain = n-kk;
 #pragma vector aligned
 #pragma ivdep
 for ( size_type k = 0 ; k < 16 ; ++k )
	  sum[k] =  0;

 if(doalpha!=-1)
 {
   const size_type iEntryBegin = m_A.graph.row_map(iRow);
   const size_type iEntryEnd   = m_A.graph.row_map(iRow+1);
   for ( size_type iEntry = iEntryBegin ; iEntry < iEntryEnd ; ++iEntry ) {
 	const scalar_type val = m_A.values(iEntry);
 	const size_type ind = m_A.graph.entries(iEntry);
     #pragma vector aligned
     #pragma ivdep
     #pragma unroll
     for ( size_type k = 0 ; k < remain ; k++ )
 	   sum[k] +=  val*m_x(ind,kk+k);
   }
 }
 else
 {
   const size_type iEntryBegin = m_A.graph.row_map(iRow);
   const size_type iEntryEnd   = m_A.graph.row_map(iRow+1);
   for ( size_type iEntry = iEntryBegin ; iEntry < iEntryEnd ; ++iEntry ) {
 	 const scalar_type val = m_A.values(iEntry);
 	 const size_type ind = m_A.graph.entries(iEntry);
     #pragma vector aligned
     #pragma ivdep
     for ( size_type k = 0 ; k < remain ; ++k )
 	  sum[k] -=  val*m_x(ind,kk+k);
   }
 }

 if(doalpha*doalpha!=1) {
   #pragma vector aligned
   #pragma ivdep
   for ( size_type k = 0 ; k < remain ; ++k )
     sum[k] *= alpha(kk+k);
 }

 if(dobeta == 0) {
   #pragma vector aligned
   #pragma ivdep
   for ( size_type k = 0 ; k < remain ; ++k )
     m_y(iRow,kk+k) = sum[k] ;
 } else if(dobeta == 1) {
   #pragma vector aligned
   #pragma ivdep
   for ( size_type k = 0 ; k < remain ; ++k )
     m_y(iRow,kk+k) += sum[k] ;
 } else if(dobeta == -1) {
   #pragma vector aligned
   #pragma ivdep
   for ( size_type k = 0 ; k < remain ; ++k )
     m_y(iRow,kk+k) = -m_y(iRow,kk+k) + sum[k] ;
 } else {
   #pragma vector aligned
   #pragma ivdep
   for ( size_type k = 0 ; k < remain ; ++k )
     m_y(iRow,kk+k) = beta(kk+k)*m_y(iRow,kk+k) + sum[k] ;
 }
  }
};

template<class RangeVector,class CrsMatrix,class DomainVector,class CoeffVector1,class CoeffVector2, int UNROLL, int doalpha = 1,int dobeta = 0>
struct MV_MultiplyFunctorUnroll
{
  typedef typename CrsMatrix::device_type                   device_type ;
  typedef typename CrsMatrix::ordinal_type                    size_type ;
  typedef typename CrsMatrix::scalar_type                    scalar_type ;
//  typedef typename KokkosArray::View<scalar_type* ,typename CrsMatrix::device_type > range_values;

  CoeffVector1 beta;
  CoeffVector2 alpha;
  CrsMatrix  m_A ;
  DomainVector  m_x ;
  RangeVector  m_y ;
  size_type n;


  //--------------------------------------------------------------------------

  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const size_type iRow ) const
  {
	scalar_type __attribute__((aligned(64))) sum[UNROLL];
    #pragma unroll
    for(size_type k = 0;k<UNROLL;k++)
      sum[k] = 0;

    if(doalpha!=-1)
    {
      const size_type iEntryBegin = m_A.graph.row_map(iRow);
      const size_type iEntryEnd   = m_A.graph.row_map(iRow+1);

      #pragma vector aligned
      #pragma ivdep
      for ( size_type iEntry = iEntryBegin ; iEntry < iEntryEnd ; ++iEntry ) {
    	  const scalar_type val = m_A.values(iEntry);
    	  const size_type ind = m_A.graph.entries(iEntry);
          #pragma unroll
          for(size_type k = 0;k<UNROLL;k++)
    	    sum[k] +=  val*m_x(ind,k);
      }
    }
    else
    {
      const size_type iEntryBegin = m_A.graph.row_map(iRow);
      const size_type iEntryEnd   = m_A.graph.row_map(iRow+1);
      #pragma vector aligned
      #pragma ivdep
      for ( size_type iEntry = iEntryBegin ; iEntry < iEntryEnd ; ++iEntry ) {
    	  const scalar_type val = m_A.values(iEntry);
    	  const size_type ind = m_A.graph.entries(iEntry);

		  #pragma unroll
          for(size_type k = 0;k<UNROLL;k++)
            sum[k] -=  val*m_x(ind,k);
      }
    }

    if(doalpha*doalpha!=1) {
      #pragma vector aligned
      #pragma ivdep
      #pragma unroll
      for ( size_type k = 0 ; k < UNROLL ; ++k )
        sum[k] *= alpha(k);
    }

    if(dobeta == 0) {
      #pragma vector aligned
      #pragma ivdep
	  #pragma unroll
      for ( size_type k = 0 ; k < UNROLL ; ++k )
        m_y(iRow,k) = sum[k] ;
      return;
    }
    else
    if(dobeta == 1) {
      #pragma vector aligned
      #pragma ivdep
      #pragma unroll
      for ( size_type k = 0 ; k < UNROLL ; ++k )
        m_y(iRow,k) += sum[k] ;
      return;
    }
    else
    if(dobeta == -1) {
      #pragma vector aligned
      #pragma ivdep
      #pragma unroll
      for ( size_type k = 0 ; k < UNROLL ; ++k )
        m_y(iRow,k) = -m_y(iRow,k) + sum[k] ;
      return;
    }
    else
    {
      #pragma vector aligned
      #pragma ivdep
      #pragma unroll
      for ( size_type k = 0 ; k < UNROLL ; ++k )
        m_y(iRow,k) = beta(k)*m_y(iRow,k) + sum[k] ;
      return;
    }
  }
};


template<class RangeVector,class CrsMatrix,class DomainVector,class CoeffVector1,class CoeffVector2,
         int doalpha, int dobeta>
void MV_Multiply( const CoeffVector1 &betav, const RangeVector & y, const CoeffVector2& alphav,
		const CrsMatrix & A , const DomainVector & x,int beta, int alpha)
{
    MV_MultiplyFunctor<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,doalpha,dobeta> op ;
    const typename CrsMatrix::ordinal_type nrow = A.numRows;
    op.m_A = A ;
    op.m_x = x ;
    op.m_y = y ;
    op.beta = betav;
    op.alpha = alphav;
    op.n=x.dimension(1);
    KokkosArray::parallel_for( nrow , op );
}

template<class RangeVector,class CrsMatrix,class DomainVector,class CoeffVector1,class CoeffVector2 >
void MV_Multiply( const CoeffVector1 &betav, const RangeVector & y, const CoeffVector2& alphav,
		const CrsMatrix & A , const DomainVector & x,int beta, int alpha)
{
  if(beta==0) {
    if(alpha==0)
      MV_MulScalar<RangeVector,RangeVector>( y, 0, y);
    else if(alpha==1)
   	  MV_Multiply<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,1,0>( betav, y, alphav,A ,  x, beta,  alpha);
    else if(alpha==-1)
   	  MV_Multiply<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,-1,0>( betav, y, alphav,A ,  x, beta,  alpha);
    else
   	  MV_Multiply<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,2,0>( betav, y, alphav,A ,  x, beta,  alpha);
  } else if(beta==1){
    if(alpha==0)
  	  return;
    else if(alpha==1)
   	  MV_Multiply<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,1,1>( betav, y, alphav,A ,  x, beta,  alpha);
    else if(alpha==-1)
   	  MV_Multiply<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,-1,1>( betav, y, alphav,A ,  x, beta,  alpha);
    else
   	  MV_Multiply<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,2,1>( betav, y, alphav,A ,  x, beta,  alpha);
  } else if(beta==-1) {
    if(alpha==0)
      MV_MulScalar<RangeVector,RangeVector>( y, -1, y);
    else if(alpha==1)
   	  MV_Multiply<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,1,-1>( betav, y, alphav,A ,  x, beta,  alpha);
    else if(alpha==-1)
   	  MV_Multiply<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,-1,-1>( betav, y, alphav,A ,  x, beta,  alpha);
    else
   	  MV_Multiply<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,2,-1>( betav, y, alphav,A ,  x, beta,  alpha);
  } else {
    if(alpha==0)
      MV_MulScalar( y, betav, y);
    else if(alpha==1)
   	  MV_Multiply<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,1,2>( betav, y, alphav,A ,  x, beta,  alpha);
    else if(alpha==-1)
   	  MV_Multiply<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,-1,2>( betav, y, alphav,A ,  x, beta,  alpha);
    else
   	  MV_Multiply<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,2,2>( betav, y, alphav,A ,  x, beta,  alpha);
  }
}

template<class RangeVector,class CrsMatrix,class DomainVector,class CoeffVector1,class CoeffVector2,
         int UNROLL, int doalpha, int dobeta>
void MV_MultiplyUnroll( const CoeffVector1 &betav, const RangeVector & y, const CoeffVector2& alphav,
		const CrsMatrix & A , const DomainVector & x,int beta, int alpha)
{
    MV_MultiplyFunctorUnroll<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,UNROLL,doalpha,dobeta> op ;
    const typename CrsMatrix::ordinal_type nrow = A.numRows;
    op.m_A = A ;
    op.m_x = x ;
    op.m_y = y ;
    op.beta = betav;
    op.alpha = alphav;
    op.n=x.dimension_1();
    KokkosArray::parallel_for( nrow , op );
}

template<class RangeVector,class CrsMatrix,class DomainVector,class CoeffVector1,class CoeffVector2,int UNROLL>
void MV_MultiplyUnroll( const CoeffVector1 &betav, const RangeVector & y, const CoeffVector2& alphav,
		const CrsMatrix & A , const DomainVector & x,int beta, int alpha)
{
  if(beta==0) {
    if(alpha==0)
      MV_MulScalar<RangeVector,RangeVector>( y, 0, y);
    else if(alpha==1)
   	  MV_MultiplyUnroll<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,UNROLL,1,0>( betav, y, alphav,A ,  x, beta,  alpha);
    else if(alpha==-1)
   	  MV_MultiplyUnroll<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,UNROLL,-1,0>( betav, y, alphav,A ,  x, beta,  alpha);
    else
   	  MV_MultiplyUnroll<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,UNROLL,2,0>( betav, y, alphav,A ,  x, beta,  alpha);
  } else if(beta==1){
    if(alpha==0)
      return;
    else if(alpha==1)
   	  MV_MultiplyUnroll<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,UNROLL,1,1>( betav, y, alphav,A ,  x, beta,  alpha);
    else if(alpha==-1)
   	  MV_MultiplyUnroll<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,UNROLL,-1,1>( betav, y, alphav,A ,  x, beta,  alpha);
    else
   	  MV_MultiplyUnroll<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,UNROLL,2,1>( betav, y, alphav,A ,  x, beta,  alpha);
  } else if(beta==-1) {
    if(alpha==0)
      MV_MulScalar<RangeVector,RangeVector>( y, -1, y);
    else if(alpha==1)
   	  MV_MultiplyUnroll<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,UNROLL,1,-1>( betav, y, alphav,A ,  x, beta,  alpha);
    else if(alpha==-1)
   	  MV_MultiplyUnroll<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,UNROLL,-1,-1>( betav, y, alphav,A ,  x, beta,  alpha);
    else
   	  MV_MultiplyUnroll<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,UNROLL,2,-1>( betav, y, alphav,A ,  x, beta,  alpha);
  } else {
    if(alpha==0)
      MV_MulScalar( y, betav, y);
    else if(alpha==1)
   	  MV_MultiplyUnroll<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,UNROLL,1,2>( betav, y, alphav,A ,  x, beta,  alpha);
    else if(alpha==-1)
   	  MV_MultiplyUnroll<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,UNROLL,-1,2>( betav, y, alphav,A ,  x, beta,  alpha);
    else
   	  MV_MultiplyUnroll<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,UNROLL,2,2>( betav, y, alphav,A ,  x, beta,  alpha);
  }
}

template<class RangeVector,class CrsMatrix,class DomainVector,class CoeffVector1,class CoeffVector2>
void MV_MultiplyUnroll( const CoeffVector1 &betav, const RangeVector & y, const CoeffVector2& alphav,
		const CrsMatrix & A , const DomainVector & x,int beta, int alpha)
{
	switch(x.dimension(1)) {
	  case 1:  MV_MultiplyUnroll<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,1 >( betav, y, alphav,A ,  x, beta,  alpha); break;
	  case 2:  MV_MultiplyUnroll<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,2 >( betav, y, alphav,A ,  x, beta,  alpha); break;
	  case 3:  MV_MultiplyUnroll<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,3 >( betav, y, alphav,A ,  x, beta,  alpha); break;
	  case 4:  MV_MultiplyUnroll<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,4 >( betav, y, alphav,A ,  x, beta,  alpha); break;
	  case 5:  MV_MultiplyUnroll<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,5 >( betav, y, alphav,A ,  x, beta,  alpha); break;
	  case 6:  MV_MultiplyUnroll<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,6 >( betav, y, alphav,A ,  x, beta,  alpha); break;
	  case 7:  MV_MultiplyUnroll<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,7 >( betav, y, alphav,A ,  x, beta,  alpha); break;
	  case 8:  MV_MultiplyUnroll<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,8 >( betav, y, alphav,A ,  x, beta,  alpha); break;
	  case 9:  MV_MultiplyUnroll<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,9 >( betav, y, alphav,A ,  x, beta,  alpha); break;
	  case 10: MV_MultiplyUnroll<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,10>( betav, y, alphav,A ,  x, beta,  alpha); break;
	  case 11: MV_MultiplyUnroll<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,11>( betav, y, alphav,A ,  x, beta,  alpha); break;
	  case 12: MV_MultiplyUnroll<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,12>( betav, y, alphav,A ,  x, beta,  alpha); break;
	  case 13: MV_MultiplyUnroll<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,13>( betav, y, alphav,A ,  x, beta,  alpha); break;
	  case 14: MV_MultiplyUnroll<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,14>( betav, y, alphav,A ,  x, beta,  alpha); break;
	  case 15: MV_MultiplyUnroll<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,15>( betav, y, alphav,A ,  x, beta,  alpha); break;
	  case 16: MV_MultiplyUnroll<RangeVector,CrsMatrix,DomainVector,CoeffVector1,CoeffVector2,16>( betav, y, alphav,A ,  x, beta,  alpha); break;
	}
}

template<class RangeVector,class CrsMatrix,class DomainVector
,typename Value1, class Layout1, class Device1, class MemoryManagement1
,typename Value2, class Layout2, class Device2, class MemoryManagement2>
void MV_Multiply(
		const KokkosArray::View<Value1,Layout1,Device1,MemoryManagement1>& betav, const RangeVector & y,
		const KokkosArray::View<Value2,Layout2,Device2,MemoryManagement2>& alphav,
		const CrsMatrix & A , const DomainVector & x)
{
  if(x.dimension(1)<=16)
    return MV_MultiplyUnroll(betav,y,alphav,A,x,2,2);
  return MV_Multiply(betav,y,alphav,A,x,2,2);
}

template<class RangeVector,class CrsMatrix,class DomainVector
,typename Value1, class Layout1, class Device1, class MemoryManagement1>
void MV_Multiply(const RangeVector & y,
		const KokkosArray::View<Value1,Layout1,Device1,MemoryManagement1>& alphav,
		const CrsMatrix & A , const DomainVector & x)
{
  if(x.dimension(1)<=16)
    return MV_MultiplyUnroll(alphav,y,alphav,A,x,0,2);
  return MV_Multiply(alphav,y,alphav,A,x,0,2);
}

template<class RangeVector,class CrsMatrix,class DomainVector>
void MV_Multiply( const RangeVector & y,
		const CrsMatrix & A , const DomainVector & x)
{
#ifdef KOKKOS_USE_CUSPARSE
    if(MV_Multiply_Try_CuSparse(0.0,y,1.0,A,x)) return;
#endif
    typedef KokkosArray::View<typename DomainVector::scalar_type*,typename DomainVector::device_type> aVector;
    aVector a;
    if(x.dimension(1)<=16)
     return MV_MultiplyUnroll(a,y,a,A,x,0,1);
    return MV_Multiply(a,y,a,A,x,0,1);
}

template<class RangeVector,class CrsMatrix,class DomainVector>
void MV_Multiply( const RangeVector & y, typename DomainVector::scalar_type s_a,
		const CrsMatrix & A , const DomainVector & x)
{
#ifdef KOKKOS_USE_CUSPARSE
    if(MV_Multiply_Try_CuSparse(0.0,y,s_a,A,x)) return;
#endif
    typedef KokkosArray::View<typename RangeVector::scalar_type*,typename RangeVector::device_type> aVector;
    aVector a;
    const int numVecs = x.dimension_1();
    if(s_a==-1) {
      if(numVecs<=UNROLL_LIMIT)
        return MV_MultiplyUnroll(a,y,a,A,x,0,-1);
      return MV_Multiply(a,y,a,A,x,0,-1);
    } else  if(s_a==1) {
      if(numVecs<=16)
        return MV_MultiplyUnroll(a,y,a,A,x,0,1);
      return MV_Multiply(a,y,a,A,x,0,1);
    } if(s_a!=0) {
  	  a = aVector("a",numVecs);
  	  typename aVector::HostMirror h_a = KokkosArray::create_mirror_view(a);
  	  for(int i=0;i<numVecs;i++)
  	    h_a(i) = s_a;
  	  KokkosArray::deep_copy(a,h_a);
      if(numVecs<=16)
        return MV_MultiplyUnroll(a,y,a,A,x,0,2);
      return MV_Multiply(a,y,a,A,x,0,2);
    }
}

//ToDo: make type of s_b and s_a independent of MV scalar_types
template<class RangeVector,class CrsMatrix,class DomainVector>
void MV_Multiply( typename RangeVector::scalar_type s_b,const RangeVector & y, typename DomainVector::scalar_type s_a,
		const CrsMatrix & A , const DomainVector & x)
{
#ifdef KOKKOS_USE_CUSPARSE
    if(MV_Multiply_Try_CuSparse(s_b,y,s_a,A,x)) return;
#endif
    typedef KokkosArray::View<typename RangeVector::scalar_type*,typename RangeVector::device_type> aVector;
    aVector a;
    aVector b;
    int numVecs = x.dimension_1();
    if(x.dimension(1)<=UNROLL_LIMIT){
      if(s_b==0) {
        if(s_a==0)
    	  return MV_MultiplyUnroll(a,y,a,A,x,0,0);
        else if(s_a==1)
      	  return MV_MultiplyUnroll(a,y,a,A,x,0,1);
        else if(s_a==-1)
      	  return MV_MultiplyUnroll(a,y,a,A,x,0,-1);
        else{
    	  a = aVector("a",numVecs);
    	  typename aVector::HostMirror h_a = KokkosArray::create_mirror_view(a);
    	  for(int i=0;i<numVecs;i++)
    	    h_a(i) = s_a;
    	  KokkosArray::deep_copy(a,h_a);
      	  return MV_MultiplyUnroll(a,y,a,A,x,0,2);
        }
      } else if(s_b==1) {
        if(s_a==0)
          return MV_MultiplyUnroll(a,y,a,A,x,1,0);
        else if(s_a==1)
          return MV_MultiplyUnroll(a,y,a,A,x,1,1);
        else if(s_a==-1)
       	  return MV_MultiplyUnroll(a,y,a,A,x,1,-1);
        else{
      	  a = aVector("a",numVecs);
      	  typename aVector::HostMirror h_a = KokkosArray::create_mirror_view(a);
      	  for(int i=0;i<numVecs;i++)
      	    h_a(i) = s_a;
          KokkosArray::deep_copy(a,h_a);
          return MV_MultiplyUnroll(a,y,a,A,x,1,2);
        }
      } else if(s_b==-1) {
        if(s_a==0)
          return MV_MultiplyUnroll(a,y,a,A,x,-1,0);
        else if(s_a==1)
          return MV_MultiplyUnroll(a,y,a,A,x,-1,1);
        else if(s_a==-1)
       	  return MV_MultiplyUnroll(a,y,a,A,x,-1,-1);
        else{
       	  a = aVector("a",numVecs);
       	  typename aVector::HostMirror h_a = KokkosArray::create_mirror_view(a);
       	  for(int i=0;i<numVecs;i++)
       	    h_a(i) = s_a;
          KokkosArray::deep_copy(a,h_a);
          return MV_MultiplyUnroll(a,y,a,A,x,-1,2);
        }
      } else {
       	b = aVector("b",numVecs);
       	typename aVector::HostMirror h_b = KokkosArray::create_mirror_view(b);
       	for(int i=0;i<numVecs;i++)
     	  h_b(i) = s_b;
        KokkosArray::deep_copy(b,h_b);
        if(s_a==0)
          return MV_MultiplyUnroll(b,y,a,A,x,2,0);
        else if(s_a==1)
          return MV_MultiplyUnroll(b,y,a,A,x,2,1);
        else if(s_a==-1)
       	  return MV_MultiplyUnroll(b,y,a,A,x,2,-1);
        else{
       	  a = aVector("a",numVecs);
       	  typename aVector::HostMirror h_a = KokkosArray::create_mirror_view(a);
       	  for(int i=0;i<numVecs;i++)
       	    h_a(i) = s_a;
          KokkosArray::deep_copy(a,h_a);
          return MV_MultiplyUnroll(b,y,a,A,x,2,2);
        }
      }
    } else {


        if(s_b==0) {
          if(s_a==0)
      	  return MV_Multiply(a,y,a,A,x,0,0);
          else if(s_a==1)
        	  return MV_Multiply(a,y,a,A,x,0,1);
          else if(s_a==-1)
        	  return MV_Multiply(a,y,a,A,x,0,-1);
          else{
      	    a = aVector("a",numVecs);
      	    typename aVector::HostMirror h_a = KokkosArray::create_mirror_view(a);
      	    for(int i=0;i<numVecs;i++)
      	      h_a(i) = s_a;
      	    KokkosArray::deep_copy(a,h_a);
        	return MV_Multiply(a,y,a,A,x,0,2);
          }
        } else if(s_b==1) {
          if(s_a==0)
            return MV_Multiply(a,y,a,A,x,1,0);
          else if(s_a==1)
            return MV_Multiply(a,y,a,A,x,1,1);
          else if(s_a==-1)
         	  return MV_Multiply(a,y,a,A,x,1,-1);
          else{
        	  a = aVector("a",numVecs);
        	  typename aVector::HostMirror h_a = KokkosArray::create_mirror_view(a);
        	  for(int i=0;i<numVecs;i++)
        	    h_a(i) = s_a;
            KokkosArray::deep_copy(a,h_a);
            return MV_Multiply(a,y,a,A,x,1,2);
          }
        } else if(s_b==-1) {
          if(s_a==0)
            return MV_Multiply(a,y,a,A,x,-1,0);
          else if(s_a==1)
            return MV_Multiply(a,y,a,A,x,-1,1);
          else if(s_a==-1)
         	  return MV_Multiply(a,y,a,A,x,-1,-1);
          else{
         	  a = aVector("a",numVecs);
         	  typename aVector::HostMirror h_a = KokkosArray::create_mirror_view(a);
         	  for(int i=0;i<numVecs;i++)
         	    h_a(i) = s_a;
            KokkosArray::deep_copy(a,h_a);
            return MV_Multiply(a,y,a,A,x,-1,2);
          }
        } else {
          b = aVector("b",numVecs);
          typename aVector::HostMirror h_b = KokkosArray::create_mirror_view(b);
          for(int i=0;i<numVecs;i++)
       	    h_b(i) = s_b;
          KokkosArray::deep_copy(b,h_b);
          if(s_a==0)
            return MV_Multiply(b,y,a,A,x,2,0);
          else if(s_a==1)
            return MV_Multiply(b,y,a,A,x,2,1);
          else if(s_a==-1)
         	  return MV_Multiply(b,y,a,A,x,2,-1);
          else{
         	  a = aVector("a",numVecs);
         	  typename aVector::HostMirror h_a = KokkosArray::create_mirror_view(a);
         	  for(int i=0;i<numVecs;i++)
         	    h_a(i) = s_a;
            KokkosArray::deep_copy(a,h_a);
            return MV_Multiply(b,y,a,A,x,2,2);
          }
        }
    }
}

namespace KokkosCrsMatrix
{
template <class CrsMatrixDst,class CrsMatrixSrc>
void deep_copy(CrsMatrixDst A,CrsMatrixSrc B)
{
  KokkosArray::deep_copy(A.graph.entries,B.graph.entries);
  //KokkosArray::deep_copy(A.graph.row_map,B.graph.row_map);
  KokkosArray::deep_copy(A.values,B.values);
}
}

}
#endif /* KOKKOSARRAY_CRSMATRIX_H_ */
