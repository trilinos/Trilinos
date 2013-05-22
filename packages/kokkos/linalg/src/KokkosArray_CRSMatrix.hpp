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
#ifdef KOKKOS_USE_MKL
#include <mkl.h>
#include <mkl_spblas.h>
#include <Kokkos_CRSMatrix_MKL.hpp>
#endif

#ifndef KOKKOSARRAY_INLINE_FUNCTION
#define KOKKOSARRAY_INLINE_FUNCTION
#endif

namespace KokkosArray
{

template<class MatrixType>
struct SparseRowView {
  typedef typename MatrixType::scalar_type scalar_type;
  typedef typename MatrixType::ordinal_type ordinal_type;

private:
  scalar_type* _values;
  ordinal_type* _colidx;
  const int _stride;

public:
  KOKKOSARRAY_INLINE_FUNCTION
  SparseRowView(scalar_type* const& values, ordinal_type* const& colidx, const int& stride, const int& count):
	  _values(values), _colidx(colidx), _stride(stride), length(count) {
  }

  const int length;

  KOKKOSARRAY_INLINE_FUNCTION
  scalar_type & value(const int& i) const {
    return _values[i*_stride];
  }

  KOKKOSARRAY_INLINE_FUNCTION
  ordinal_type & colidx(const int& i) const {
    return _colidx[i*_stride];
  }


};

template< typename ScalarType , typename OrdinalType, class Device, class MemoryTraits = void >
class CrsMatrix
{
  public:
    typedef Device      device_type ;
    typedef ScalarType  scalar_type ;
    typedef OrdinalType  ordinal_type ;

#ifdef _OPENMP
    typedef KokkosArray::OpenMP host_device_type;
#else
    typedef KokkosArray::Host host_device_type;
#endif
    typedef CrsMatrix<ScalarType, OrdinalType, host_device_type, MemoryTraits> HostMirror;
    typedef KokkosArray::CrsArray<OrdinalType, KokkosArray::LayoutLeft, Device,int> CrsArrayType;
    typedef typename CrsArrayType::entries_type   index_type;
    typedef typename CrsArrayType::row_map_type   row_map_type;
    typedef KokkosArray::View< scalar_type* , KokkosArray::LayoutLeft, device_type,MemoryTraits >   values_type ;

    typedef typename values_type::const_scalar_type  const_scalar_type;
    typedef typename values_type::non_const_scalar_type  non_const_scalar_type;


#ifdef KOKKOS_USE_CUSPARSE
    cusparseHandle_t cusparse_handle;
    cusparseMatDescr_t cusparse_descr;
#endif
    CrsArrayType graph;
    values_type  values;

    CrsMatrix() {
      _numRows = 0;
      _numCols = 0;
      _nnz = 0;
    };

    CrsMatrix(const std::string &label, OrdinalType nrows, OrdinalType ncols, OrdinalType annz, ScalarType* val, OrdinalType* rows, OrdinalType* cols, bool pad = false) {
      import(label, nrows, ncols, annz, val, rows, cols);
#ifdef KOKKOS_USE_CUSPARSE
      cusparseCreate(&cusparse_handle);
      cusparseCreateMatDescr(&cusparse_descr);
#endif
    };

    CrsMatrix(const std::string &label, OrdinalType nrows, OrdinalType ncols, OrdinalType annz, values_type vals, row_map_type rows, index_type cols) {
    	_numRows = nrows;
    	_numCols = ncols;
    	_nnz = annz;
    	graph.row_map = rows;
    	graph.entries = cols;
    	values = vals;
#ifdef KOKKOS_USE_CUSPARSE
      cusparseCreate(&cusparse_handle);
      cusparseCreateMatDescr(&cusparse_descr);
#endif
    };

    void import(const std::string &label, OrdinalType nrows, OrdinalType ncols, OrdinalType annz, ScalarType* val, OrdinalType* rows, OrdinalType* cols);
    void generate(const std::string &label, OrdinalType nrows, OrdinalType ncols, OrdinalType target_nnz, OrdinalType varianz_nel_row, OrdinalType width_row);

    template< typename aScalarType , typename aOrdinalType, class aDevice, class aMemoryTraits>
    CrsMatrix& operator = (const CrsMatrix<aScalarType,aOrdinalType,aDevice,aMemoryTraits> & mtx) {
      _numRows = mtx.numRows();
      _numCols = mtx.numCols();
      _nnz = mtx.nnz();
      graph = mtx.graph;
      values = mtx.values;
      return *this;
    }

    KOKKOSARRAY_INLINE_FUNCTION
    ordinal_type numRows() const
       {return _numRows;};
    KOKKOSARRAY_INLINE_FUNCTION
    ordinal_type numCols() const
       {return _numCols;};
    KOKKOSARRAY_INLINE_FUNCTION
    ordinal_type nnz() const
       {return _nnz;};

    friend struct SparseRowView<CrsMatrix>;
    KOKKOSARRAY_INLINE_FUNCTION
    SparseRowView<CrsMatrix> row(int i) const{
      const int start = graph.row_map(i);
      const int end = graph.row_map(i+1);
      return SparseRowView<CrsMatrix>(&values(start),&graph.entries(start),1,end-start);
    }
  private:
    ordinal_type _numRows;
    ordinal_type _numCols;
    ordinal_type _nnz;


};

template< typename ScalarType , typename OrdinalType, class Device, class MemoryTraits >
void CrsMatrix<ScalarType , OrdinalType, Device, MemoryTraits >::import(const std::string &label, OrdinalType nrows, OrdinalType ncols, OrdinalType annz, ScalarType* val, OrdinalType* rows, OrdinalType* cols)
{
  std::string str = label;
  values = values_type(str.append(".values"), annz);

  _numRows = nrows;
  _numCols = ncols;
  _nnz = annz;

  std::vector<int> row_lengths(_numRows, 0);

  for(int i = 0; i < _numRows ; i++)
    row_lengths[i] = rows[i + 1] - rows[i];

  str = label;
  graph = KokkosArray::create_crsarray<CrsArrayType>(str.append(".graph"), row_lengths);
  typename values_type::HostMirror h_values = KokkosArray::create_mirror_view(values);
  typename index_type::HostMirror h_entries = KokkosArray::create_mirror_view(graph.entries);

  for(OrdinalType i = 0; i < _nnz; i++) {
    h_values(i) = val[i];
    h_entries(i) = cols[i];
  }

  KokkosArray::deep_copy(values, h_values);
  KokkosArray::deep_copy(graph.entries, h_entries);
}

template< typename ScalarType , typename OrdinalType, class Device, class MemoryTraits >
void CrsMatrix<ScalarType , OrdinalType, Device, MemoryTraits >::generate(const std::string &label, OrdinalType nrows, OrdinalType ncols, OrdinalType target_nnz, OrdinalType varianz_nel_row, OrdinalType width_row)
{
  _numRows = nrows;
  _numCols = ncols;

  graph.row_map = row_map_type("CrsMatrix::rowPtr", nrows + 1);
  typename row_map_type::HostMirror h_row_map = KokkosArray::create_mirror_view(graph.row_map);

  OrdinalType elements_per_row = target_nnz / nrows;
  srand(13721);
  h_row_map(0) = 0;

  for(int row = 0; row < nrows; row++) {
   // int varianz = (1.0 * rand() / INT_MAX - 0.5) * varianz_nel_row;
   // h_row_map(row + 1) = h_row_map(row) + elements_per_row + varianz;
  }

  _nnz = h_row_map(nrows);
  values = values_type("CrsMatrix::values", _nnz);
  graph.entries = index_type("CrsMatrix::colInd", _nnz);
  typename values_type::HostMirror h_values = KokkosArray::create_mirror_view(values);
  typename index_type::HostMirror h_entries = KokkosArray::create_mirror_view(graph.entries);

  for(int row = 0; row < nrows; row++) {
    for(int k = h_row_map(row); k < h_row_map(row + 1); row++) {
      //int pos = (1.0 * rand() / INT_MAX - 0.5) * width_row;

      //if(pos < 0) pos += ncols;

     // if(pos >= ncols) pos -= ncols;

     // h_entries(k) = pos;
     // h_values(k) = 100.0 * rand() / INT_MAX - 50.0;
    }
  }

  KokkosArray::deep_copy(values, h_values);
  KokkosArray::deep_copy(graph.entries, h_entries);
  KokkosArray::deep_copy(graph.row_map, h_row_map);
}

template<typename Scalar>
KOKKOSARRAY_INLINE_FUNCTION
Scalar shfl_down(const Scalar &val, const int& delta, const int& width){
  return val;
}

template<>
KOKKOSARRAY_INLINE_FUNCTION
unsigned int shfl_down<unsigned int>(const unsigned int &val, const int& delta, const int& width){
#ifdef __CUDA_ARCH__
  #if (__CUDA_ARCH__ >= 300)
    unsigned int tmp1 = val;
    int tmp = *reinterpret_cast<int*>(&tmp1);
    tmp = __shfl_down(tmp,delta,width);
    return *reinterpret_cast<unsigned int*>(&tmp);
  #else
    return val;
  #endif
#else
  return val;
#endif
}

template<>
KOKKOSARRAY_INLINE_FUNCTION
int shfl_down<int>(const int &val, const int& delta, const int& width){
#ifdef __CUDA_ARCH__
  #if (__CUDA_ARCH__ >= 300)
    return __shfl_down(val,delta,width);
  #else
    return val;
  #endif
#else
  return val;
#endif
}

template<>
KOKKOSARRAY_INLINE_FUNCTION
float shfl_down<float>(const float &val, const int& delta, const int& width){
#ifdef __CUDA_ARCH__
  #if (__CUDA_ARCH__ >= 300)
    return __shfl_down(val,delta,width);
  #else
    return val;
  #endif
#else
  return val;
#endif
}

template<>
KOKKOSARRAY_INLINE_FUNCTION
double shfl_down<double>(const double &val, const int& delta, const int& width){
#ifdef __CUDA_ARCH__
  #if (__CUDA_ARCH__ >= 300)
    int lo = __double2loint(val);
    int hi = __double2hiint(val);
    lo = __shfl_down(lo,delta,width);
    hi = __shfl_down(hi,delta,width);
    return __hiloint2double(hi,lo);
  #else
    return val;
  #endif
#else
  return val;
#endif
}


template<class RangeVector, class CrsMatrix, class DomainVector, class CoeffVector1, class CoeffVector2, int doalpha , int dobeta, int ThreadsPerRow>
struct MV_MultiplyFunctor {
  typedef typename CrsMatrix::device_type                   device_type ;
  typedef typename CrsMatrix::ordinal_type                    size_type ;
  typedef typename CrsMatrix::non_const_scalar_type         scalar_type ;
  typedef typename KokkosArray::View<scalar_type* , typename CrsMatrix::device_type > range_values;

  CoeffVector1 beta;
  CoeffVector2 alpha;
  CrsMatrix  m_A ;
  DomainVector  m_x ;
  RangeVector  m_y ;
  size_type n;

  //--------------------------------------------------------------------------


  template<int UNROLL>
  KOKKOSARRAY_INLINE_FUNCTION
  void strip_mine(const size_type i, const size_type kk) const {
    const size_type iRow = i/ThreadsPerRow;
    const int lane = i%ThreadsPerRow;
    scalar_type sum[UNROLL];
#pragma ivdep
#pragma unroll
    for(size_type k = 0 ; k < UNROLL ; ++k)
      sum[k] =  0;

    if(doalpha != -1) {
  	  const SparseRowView<CrsMatrix> row = m_A.row(iRow);

      #pragma ivdep
  	  #pragma loop count (24)
	  #pragma unroll
      for(size_type iEntry = lane ; iEntry < row.length ; iEntry+=ThreadsPerRow) {
        const scalar_type val = row.value(iEntry);
        const size_type ind = row.colidx(iEntry);

        #pragma unroll
        for(size_type k = 0 ; k < UNROLL ; k++)
          sum[k] +=  val * m_x(ind, kk + k);
      }
    } else {
      const SparseRowView<CrsMatrix> row = m_A.row(iRow);

      #pragma ivdep
      #pragma loop count (24)
	  #pragma unroll
      for(size_type iEntry = lane ; iEntry < row.length ; iEntry+=ThreadsPerRow) {
        const scalar_type val = row.value(iEntry);
        const size_type ind = row.colidx(iEntry);

        #pragma unroll
        for(size_type k = 0 ; k < UNROLL ; ++k)
          sum[k] -=  val * m_x(ind, kk + k);
      }
    }
    for(int i=0;i<UNROLL;i++) {
if(ThreadsPerRow > 1)
    sum[i] += shfl_down(sum[i], 1,ThreadsPerRow);
if(ThreadsPerRow > 2)
	sum[i] += shfl_down(sum[i], 2,ThreadsPerRow);
if(ThreadsPerRow > 4)
	sum[i] += shfl_down(sum[i], 4,ThreadsPerRow);
if(ThreadsPerRow > 8)
	sum[i] += shfl_down(sum[i], 8,ThreadsPerRow);
if(ThreadsPerRow > 16)
	sum[i] += shfl_down(sum[i], 16,ThreadsPerRow);
    }
	if(lane==0) {

    if(doalpha * doalpha != 1) {
#pragma ivdep
#pragma unroll
      for(size_type k = 0 ; k < UNROLL ; ++k)
        sum[k] *= alpha(kk + k);
    }

    if(dobeta == 0) {
#pragma ivdep
#pragma unroll
      for(size_type k = 0 ; k < UNROLL ; ++k)
        m_y(iRow, kk + k) = sum[k] ;
    } else if(dobeta == 1) {
#pragma ivdep
#pragma unroll
      for(size_type k = 0 ; k < UNROLL ; ++k)
        m_y(iRow, kk + k) += sum[k] ;
    } else if(dobeta == -1) {
#pragma ivdep
#pragma unroll
      for(size_type k = 0 ; k < UNROLL ; ++k)
        m_y(iRow, kk + k) = -m_y(iRow, kk + k) +  sum[k] ;
    } else {
#pragma ivdep
#pragma unroll
      for(size_type k = 0 ; k < UNROLL ; ++k)
        m_y(iRow, kk + k) = beta(kk + k) * m_y(iRow, kk + k) + sum[k] ;
    }
	}
  }

  KOKKOSARRAY_INLINE_FUNCTION
  void strip_mine_1(const size_type i) const {
    const size_type iRow = i/ThreadsPerRow;
    const int lane = i%ThreadsPerRow;
    scalar_type sum = 0;

    if(doalpha != -1) {
	  const SparseRowView<CrsMatrix> row = m_A.row(iRow);

	  #pragma loop count (24)
      for(size_type iEntry = lane ; iEntry < row.length ; iEntry+=ThreadsPerRow) {
        sum +=  row.value(iEntry) * m_x(row.colidx(iEntry),0);
      }
    } else {
	  const SparseRowView<CrsMatrix> row = m_A.row(iRow);

	  #pragma loop count (24)
      for(size_type iEntry = lane ; iEntry < row.length ; iEntry+=ThreadsPerRow) {
        sum -=  row.value(iEntry) * m_x(row.colidx(iEntry),0);
      }
    }

    if(ThreadsPerRow > 1)
        sum += shfl_down(sum, 1,ThreadsPerRow);
    if(ThreadsPerRow > 2)
    	sum += shfl_down(sum, 2,ThreadsPerRow);
    if(ThreadsPerRow > 4)
    	sum += shfl_down(sum, 4,ThreadsPerRow);
    if(ThreadsPerRow > 8)
    	sum += shfl_down(sum, 8,ThreadsPerRow);
    if(ThreadsPerRow > 16)
    	sum += shfl_down(sum, 16,ThreadsPerRow);

	if(lane==0) {
    if(doalpha * doalpha != 1) {
      sum *= alpha(0);
    }

    if(dobeta == 0) {
      m_y(iRow, 0) = sum ;
    } else if(dobeta == 1) {
        m_y(iRow, 0) += sum ;
    } else if(dobeta == -1) {
        m_y(iRow, 0) = -m_y(iRow, 0) +  sum;
    } else {
      m_y(iRow, 0) = beta(0) * m_y(iRow, 0) + sum;
    }
	}
  }


  KOKKOSARRAY_INLINE_FUNCTION
  void operator()(const size_type iRow) const {
    size_type kk = 0;

#ifdef KOKKOS_FAST_COMPILE
    for(; kk + 4 <= n; kk += 4)
      strip_mine<4>(iRow, kk);
    for( ; kk<n; kk++)
      strip_mine<1>(iRow, kk);
#else
#if (DEVICE == 2)
    if((n > 8) && (n % 8 == 1)) {
      strip_mine<9>(iRow, kk);
      kk += 9;
    }
    for(; kk + 8 <= n; kk += 8)
      strip_mine<8>(iRow, kk);
    if(kk < n)
      switch(n - kk) {
#else
    if((n > 16) && (n % 16 == 1)) {
      strip_mine<17>(iRow, kk);
      kk += 17;
    }

    for(; kk + 16 <= n; kk += 16)
      strip_mine<16>(iRow, kk);

    if(kk < n)
      switch(n - kk) {
        case 15:
          strip_mine<15>(iRow, kk);
          break;

        case 14:
          strip_mine<14>(iRow, kk);
          break;

        case 13:
          strip_mine<13>(iRow, kk);
          break;

        case 12:
          strip_mine<12>(iRow, kk);
          break;

        case 11:
          strip_mine<11>(iRow, kk);
          break;

        case 10:
          strip_mine<10>(iRow, kk);
          break;

        case 9:
          strip_mine<9>(iRow, kk);
          break;

        case 8:
          strip_mine<8>(iRow, kk);
          break;
#endif
        case 7:
          strip_mine<7>(iRow, kk);
          break;

        case 6:
          strip_mine<6>(iRow, kk);
          break;

        case 5:
          strip_mine<5>(iRow, kk);
          break;

        case 4:
          strip_mine<4>(iRow, kk);
          break;

        case 3:
          strip_mine<3>(iRow, kk);
          break;

        case 2:
          strip_mine<2>(iRow, kk);
          break;

        case 1:
          strip_mine_1(iRow);
          break;
      }
#endif
  }
};

template<class RangeVector, class CrsMatrix, class DomainVector, class CoeffVector1, class CoeffVector2, int doalpha, int dobeta, int ThreadsPerRow>
struct MV_MultiplySingleFunctor {
  typedef typename CrsMatrix::device_type                   device_type ;
  typedef typename CrsMatrix::ordinal_type                    size_type ;
  typedef typename CrsMatrix::non_const_scalar_type         scalar_type ;
  typedef typename KokkosArray::View<scalar_type* , typename CrsMatrix::device_type > range_values;

  CoeffVector1 beta;
  CoeffVector2 alpha;
  CrsMatrix  m_A ;
  DomainVector  m_x ;
  RangeVector  m_y ;
  size_type n;

  //--------------------------------------------------------------------------



  KOKKOSARRAY_INLINE_FUNCTION
  void operator()(const size_type i) const {
	    const size_type iRow = i/ThreadsPerRow;
	    const int lane = i%ThreadsPerRow;
	    scalar_type sum = 0;


		if(doalpha != -1) {
		  const SparseRowView<CrsMatrix> row = m_A.row(iRow);

		  #pragma loop count (24)
		  #pragma unroll
	      for(size_type iEntry = lane ; iEntry < row.length ; iEntry+=ThreadsPerRow) {
	        sum +=  row.value(iEntry) * m_x(row.colidx(iEntry));
	      }
	    } else {
		  const SparseRowView<CrsMatrix> row = m_A.row(iRow);

          #pragma loop count (24)
          #pragma unroll
	      for(size_type iEntry = lane ; iEntry < row.length ; iEntry+=ThreadsPerRow) {
		    sum -=  row.value(iEntry) * m_x(row.colidx(iEntry));
		  }
	    }
	    if(ThreadsPerRow > 1)
	        sum += shfl_down(sum, 1,ThreadsPerRow);
	    if(ThreadsPerRow > 2)
	    	sum += shfl_down(sum, 2,ThreadsPerRow);
	    if(ThreadsPerRow > 4)
	    	sum += shfl_down(sum, 4,ThreadsPerRow);
	    if(ThreadsPerRow > 8)
	    	sum += shfl_down(sum, 8,ThreadsPerRow);
	    if(ThreadsPerRow > 16)
	    	sum += shfl_down(sum, 16,ThreadsPerRow);

		if(lane==0) {

	    if(doalpha * doalpha != 1) {
	      sum *= alpha(0);
	    }

	    if(dobeta == 0) {
	      m_y(iRow) = sum ;
	    } else if(dobeta == 1) {
	        m_y(iRow) += sum ;
	    } else if(dobeta == -1) {
	        m_y(iRow) = -m_y(iRow) +  sum;
	    } else {
	      m_y(iRow) = beta(0) * m_y(iRow) + sum;
	    }
		}
  }
};

namespace Impl{

template < class RangeVector, class CrsMatrix, class DomainVector, class CoeffVector1, class CoeffVector2>
void MV_Multiply_Check_Compatibility(const CoeffVector1 &betav, const RangeVector &y, const CoeffVector2 &alphav,
        const CrsMatrix &A , const DomainVector &x, const int & doalpha, const int & dobeta){
  typename DomainVector::size_type numVecs = x.dimension_1();
  typename DomainVector::size_type numRows = A.numRows();
  typename DomainVector::size_type numCols = A.numCols();

  if(y.dimension_1()!=numVecs) {
	  std::ostringstream msg;
	  msg << "Error in CRSMatrix - Vector Multiply (y = by + aAx): 2nd dimensions of y and x do not match\n";
	  msg << "\t Labels are: y(" << RangeVector::memory_space::query_label(y.ptr_on_device()) << ") b("
			                     << CoeffVector1::memory_space::query_label(betav.ptr_on_device()) << ") a("
			                     << CoeffVector2::memory_space::query_label(alphav.ptr_on_device()) << ") x("
			                     << CrsMatrix::values_type::memory_space::query_label(A.values.ptr_on_device()) << ") x("
			                     << DomainVector::memory_space::query_label(x.ptr_on_device()) << ")\n";
	  msg << "\t Dimensions are: y(" << y.dimension_0() << "," << y.dimension_1() << ") x(" << x.dimension_0() << "," << x.dimension_1() << ")\n";
	  Impl::throw_runtime_exception( msg.str() );
  }
  if(numRows>y.dimension_0()) {
	  std::ostringstream msg;
	  msg << "Error in CRSMatrix - Vector Multiply (y = by + aAx): dimensions of y and A do not match\n";
	  msg << "\t Labels are: y(" << RangeVector::memory_space::query_label(y.ptr_on_device()) << ") b("
			                     << CoeffVector1::memory_space::query_label(betav.ptr_on_device()) << ") a("
			                     << CoeffVector2::memory_space::query_label(alphav.ptr_on_device()) << ") x("
			                     << CrsMatrix::values_type::memory_space::query_label(A.values.ptr_on_device()) << ") x("
			                     << DomainVector::memory_space::query_label(x.ptr_on_device()) << ")\n";
	  msg << "\t Dimensions are: y(" << y.dimension_0() << "," << y.dimension_1() << ") A(" << A.numCols() << "," << A.numRows() << ")\n";
	  Impl::throw_runtime_exception( msg.str() );
  }
  if(numCols>x.dimension_0()) {
	  std::ostringstream msg;
	  msg << "Error in CRSMatrix - Vector Multiply (y = by + aAx): dimensions of x and A do not match\n";
	  msg << "\t Labels are: y(" << RangeVector::memory_space::query_label(y.ptr_on_device()) << ") b("
			                     << CoeffVector1::memory_space::query_label(betav.ptr_on_device()) << ") a("
			                     << CoeffVector2::memory_space::query_label(alphav.ptr_on_device()) << ") x("
			                     << CrsMatrix::values_type::memory_space::query_label(A.values.ptr_on_device()) << ") x("
			                     << DomainVector::memory_space::query_label(x.ptr_on_device()) << ")\n";
	  msg << "\t Dimensions are: x(" << x.dimension_0() << "," << x.dimension_1() << ") A(" << A.numCols() << "," << A.numRows() << ")\n";
	  Impl::throw_runtime_exception( msg.str() );
  }
  if(dobeta==2)
  if(betav.dimension_0()!=numVecs) {
	  std::ostringstream msg;
	  msg << "Error in CRSMatrix - Vector Multiply (y = by + aAx): 2nd dimensions of y and b do not match\n";
	  msg << "\t Labels are: y(" << RangeVector::memory_space::query_label(y.ptr_on_device()) << ") b("
			                     << CoeffVector1::memory_space::query_label(betav.ptr_on_device()) << ") a("
			                     << CoeffVector2::memory_space::query_label(alphav.ptr_on_device()) << ") x("
			                     << CrsMatrix::values_type::memory_space::query_label(A.values.ptr_on_device()) << ") x("
			                     << DomainVector::memory_space::query_label(x.ptr_on_device()) << ")\n";
	  msg << "\t Dimensions are: y(" << y.dimension_0() << "," << y.dimension_1() << ") b(" << betav.dimension_0() << ")\n";
	  Impl::throw_runtime_exception( msg.str() );
  }
  if(doalpha==2)
  if(alphav.dimension_0()!=numVecs) {
	  std::ostringstream msg;
	  msg << "Error in CRSMatrix - Vector Multiply (y = by + aAx): 2nd dimensions of x and b do not match\n";
	  msg << "\t Labels are: y(" << RangeVector::memory_space::query_label(y.ptr_on_device()) << ") b("
			                     << CoeffVector1::memory_space::query_label(betav.ptr_on_device()) << ") a("
			                     << CoeffVector2::memory_space::query_label(alphav.ptr_on_device()) << ") x("
			                     << CrsMatrix::values_type::memory_space::query_label(A.values.ptr_on_device()) << ") x("
			                     << DomainVector::memory_space::query_label(x.ptr_on_device()) << ")\n";
	  msg << "\t Dimensions are: x(" << x.dimension_0() << "," << x.dimension_1() << ") b(" << betav.dimension_0() << ")\n";
	  Impl::throw_runtime_exception( msg.str() );
  }
}
}

template<class Device, class ScalarType>
struct ThreadsPerRow {
  static const int value = 1;
};

template<>
struct ThreadsPerRow<KokkosArray::Cuda,unsigned int> {
  static const int value = 8;
};

template<>
struct ThreadsPerRow<KokkosArray::Cuda,int> {
  static const int value = 8;
};

template<>
struct ThreadsPerRow<KokkosArray::Cuda,float> {
  static const int value = 8;
};

template<>
struct ThreadsPerRow<KokkosArray::Cuda,double> {
  static const int value = 8;
};

template < class RangeVector, class TCrsMatrix, class DomainVector, class CoeffVector1, class CoeffVector2,
         int doalpha, int dobeta >
void MV_Multiply(typename KokkosArray::Impl::enable_if< DomainVector::Rank == 2, const CoeffVector1 >::type &betav, const RangeVector &y, const CoeffVector2 &alphav,
                 const TCrsMatrix &A , const DomainVector &x)
{
  if(doalpha==0) {
	  if(dobeta==2) MV_MulScalar(y,betav,y);
	  else MV_MulScalar(y,dobeta,y);
	  return;
  }
  else {
  typedef View< typename RangeVector::non_const_data_type ,
	            typename RangeVector::array_layout ,
	            typename RangeVector::device_type ,
	            typename RangeVector::memory_traits >
  RangeVectorType;

  typedef View< typename DomainVector::const_data_type ,
	            typename DomainVector::array_layout ,
	            typename DomainVector::device_type ,
	            KokkosArray::MemoryRandomRead >
  DomainVectorType;

  typedef View< typename CoeffVector1::const_data_type ,
	            typename CoeffVector1::array_layout ,
	            typename CoeffVector1::device_type ,
	            KokkosArray::MemoryRandomRead >
  CoeffVector1Type;

  typedef View< typename CoeffVector2::const_data_type ,
	            typename CoeffVector2::array_layout ,
	            typename CoeffVector2::device_type ,
	            KokkosArray::MemoryRandomRead >
  CoeffVector2Type;

  typedef CrsMatrix<typename TCrsMatrix::const_scalar_type,
		            typename TCrsMatrix::ordinal_type,
		            typename TCrsMatrix::device_type> CrsMatrixType;

  Impl::MV_Multiply_Check_Compatibility(betav,y,alphav,A,x,doalpha,dobeta);
#ifndef KOKKOS_FAST_COMPILE
  MV_MultiplyFunctor<RangeVectorType, CrsMatrixType, DomainVectorType, CoeffVector1Type, CoeffVector2Type,
  doalpha, dobeta,ThreadsPerRow<typename TCrsMatrix::device_type,typename TCrsMatrix::non_const_scalar_type>::value> op ;
  const typename CrsMatrixType::ordinal_type nrow = A.numRows();
  op.m_A = A ;
  op.m_x = x ;
  op.m_y = y ;
  op.beta = betav;
  op.alpha = alphav;
  op.n = x.dimension(1);
  KokkosArray::parallel_for(nrow*ThreadsPerRow<typename TCrsMatrix::device_type,typename TCrsMatrix::non_const_scalar_type>::value , op);

#else

  MV_MultiplyFunctor<RangeVectorType, CrsMatrixType, DomainVectorType, CoeffVector1Type, CoeffVector2Type, 2,2,
  ThreadsPerRow<typename TCrsMatrix::device_type,typename TCrsMatrix::non_const_scalar_type>::value> op ;

  int numVecs = x.dimension_1();
  CoeffVector1 beta = betav;
  CoeffVector2 alpha = alphav;

  if(doalpha!=2) {
    alpha = CoeffVector2("CrsMatrix::auto_a", numVecs);
    typename CoeffVector2::HostMirror h_a = KokkosArray::create_mirror_view(alpha);
    typename CoeffVector2::scalar_type s_a = (typename CoeffVector2::scalar_type) doalpha;
    for(int i = 0; i < numVecs; i++)
      h_a(i) = s_a;
    KokkosArray::deep_copy(alpha, h_a);
  }
  if(dobeta!=2) {
    beta = CoeffVector1("CrsMatrix::auto_b", numVecs);
    typename CoeffVector1::HostMirror h_b = KokkosArray::create_mirror_view(beta);
    typename CoeffVector1::scalar_type s_b = (typename CoeffVector1::scalar_type) dobeta;
    for(int i = 0; i < numVecs; i++)
      h_b(i) = s_b;
    KokkosArray::deep_copy(beta, h_b);
  }
  const typename CrsMatrixType::ordinal_type nrow = A.numRows();
  op.m_A = A ;
  op.m_x = x ;
  op.m_y = y ;
  op.beta = beta;
  op.alpha = alpha;
  op.n = x.dimension_1();
  KokkosArray::parallel_for(nrow*ThreadsPerRow<typename TCrsMatrix::device_type,typename TCrsMatrix::non_const_scalar_type>::value , op);
#endif
  }
}

template < class RangeVector, class TCrsMatrix, class DomainVector, class CoeffVector1, class CoeffVector2,
         int doalpha, int dobeta >
void MV_Multiply(typename KokkosArray::Impl::enable_if< DomainVector::Rank == 1, const CoeffVector1 >::type &betav, const RangeVector &y, const CoeffVector2 &alphav,
                 const TCrsMatrix &A , const DomainVector &x)
{
  if(doalpha==0) {
	  if(dobeta==2) V_MulScalar(y,betav,y);
	  else V_MulScalar(y,dobeta,y);
	  return;
  }
  else {
  typedef View< typename RangeVector::non_const_data_type ,
	            typename RangeVector::array_layout ,
	            typename RangeVector::device_type ,
	            typename RangeVector::memory_traits >
  RangeVectorType;

  typedef View< typename DomainVector::const_data_type ,
	            typename DomainVector::array_layout ,
	            typename DomainVector::device_type ,
	            KokkosArray::MemoryRandomRead >
  DomainVectorType;

  typedef View< typename CoeffVector1::const_data_type ,
	            typename CoeffVector1::array_layout ,
	            typename CoeffVector1::device_type ,
	            KokkosArray::MemoryRandomRead >
  CoeffVector1Type;

  typedef View< typename CoeffVector2::const_data_type ,
	            typename CoeffVector2::array_layout ,
	            typename CoeffVector2::device_type ,
	            KokkosArray::MemoryRandomRead >
  CoeffVector2Type;

  typedef CrsMatrix<typename TCrsMatrix::const_scalar_type,
		            typename TCrsMatrix::ordinal_type,
		            typename TCrsMatrix::device_type> CrsMatrixType;

  Impl::MV_Multiply_Check_Compatibility(betav,y,alphav,A,x,doalpha,dobeta);
#ifndef KOKKOS_FAST_COMPILE
  MV_MultiplySingleFunctor<RangeVectorType, CrsMatrixType, DomainVectorType, CoeffVector1Type, CoeffVector2Type, doalpha, dobeta
  ,ThreadsPerRow<typename CrsMatrixType::device_type,typename CrsMatrixType::non_const_scalar_type>::value> op ;
  const typename CrsMatrixType::ordinal_type nrow = A.numRows();
  op.m_A = A ;
  op.m_x = x ;
  op.m_y = y ;
  op.beta = betav;
  op.alpha = alphav;
  op.n = x.dimension(1);
  KokkosArray::parallel_for(nrow*ThreadsPerRow<typename TCrsMatrix::device_type,typename TCrsMatrix::non_const_scalar_type>::value , op);

#else

  MV_MultiplySingleFunctor<RangeVectorType, CrsMatrixType, DomainVectorType, CoeffVector1Type, CoeffVector2Type, 2,2
  ,ThreadsPerRow<typename CrsMatrixType::device_type,typename CrsMatrixType::non_const_scalar_type>::value> op ;

  int numVecs = x.dimension_1();
  CoeffVector1 beta = betav;
  CoeffVector2 alpha = alphav;

  if(doalpha!=2) {
    alpha = CoeffVector2("CrsMatrix::auto_a", numVecs);
    typename CoeffVector2::HostMirror h_a = KokkosArray::create_mirror_view(alpha);
    typename CoeffVector2::scalar_type s_a = (typename CoeffVector2::scalar_type) doalpha;
    for(int i = 0; i < numVecs; i++)
      h_a(i) = s_a;
    KokkosArray::deep_copy(alpha, h_a);
  }
  if(dobeta!=2) {
    beta = CoeffVector1("CrsMatrix::auto_b", numVecs);
    typename CoeffVector1::HostMirror h_b = KokkosArray::create_mirror_view(beta);
    typename CoeffVector1::scalar_type s_b = (typename CoeffVector1::scalar_type) dobeta;
    for(int i = 0; i < numVecs; i++)
      h_b(i) = s_b;
    KokkosArray::deep_copy(beta, h_b);
  }
  const typename CrsMatrixType::ordinal_type nrow = A.numRows();
  op.m_A = A ;
  op.m_x = x ;
  op.m_y = y ;
  op.beta = beta;
  op.alpha = alpha;
  op.n = x.dimension_1();
  KokkosArray::parallel_for(nrow*ThreadsPerRow<typename TCrsMatrix::device_type,typename TCrsMatrix::non_const_scalar_type>::value , op);
#endif
  }
}

template<class RangeVector, class CrsMatrix, class DomainVector, class CoeffVector1, class CoeffVector2 >
void MV_Multiply(const CoeffVector1 &betav, const RangeVector &y, const CoeffVector2 &alphav,
                 const CrsMatrix &A , const DomainVector &x, int beta, int alpha)
{
  if(beta == 0) {
    if(alpha == 0)
      MV_Multiply<RangeVector, CrsMatrix, DomainVector, CoeffVector1, CoeffVector2, 0, 0>(betav, y, alphav, A ,  x);
    else if(alpha == 1)
      MV_Multiply<RangeVector, CrsMatrix, DomainVector, CoeffVector1, CoeffVector2, 1, 0>(betav, y, alphav, A ,  x);
    else if(alpha == -1)
      MV_Multiply < RangeVector, CrsMatrix, DomainVector, CoeffVector1, CoeffVector2, -1, 0 > (betav, y, alphav, A ,  x);
    else
      MV_Multiply<RangeVector, CrsMatrix, DomainVector, CoeffVector1, CoeffVector2, 2, 0>(betav, y, alphav, A ,  x);
  } else if(beta == 1) {
    if(alpha == 0)
      return;
    else if(alpha == 1)
      MV_Multiply<RangeVector, CrsMatrix, DomainVector, CoeffVector1, CoeffVector2, 1, 1>(betav, y, alphav, A ,  x);
    else if(alpha == -1)
      MV_Multiply < RangeVector, CrsMatrix, DomainVector, CoeffVector1, CoeffVector2, -1, 1 > (betav, y, alphav, A ,  x);
    else
      MV_Multiply<RangeVector, CrsMatrix, DomainVector, CoeffVector1, CoeffVector2, 2, 1>(betav, y, alphav, A ,  x);
  } else if(beta == -1) {
    if(alpha == 0)
      MV_Multiply<RangeVector, CrsMatrix, DomainVector, CoeffVector1, CoeffVector2, 0, -1>(betav, y, alphav, A ,  x);
    else if(alpha == 1)
      MV_Multiply < RangeVector, CrsMatrix, DomainVector, CoeffVector1, CoeffVector2, 1, -1 > (betav, y, alphav, A ,  x);
    else if(alpha == -1)
      MV_Multiply < RangeVector, CrsMatrix, DomainVector, CoeffVector1, CoeffVector2, -1, -1 > (betav, y, alphav, A ,  x);
    else
      MV_Multiply < RangeVector, CrsMatrix, DomainVector, CoeffVector1, CoeffVector2, 2, -1 > (betav, y, alphav, A ,  x);
  } else {
    if(alpha == 0)
      MV_Multiply<RangeVector, CrsMatrix, DomainVector, CoeffVector1, CoeffVector2, 0, 2>(betav, y, alphav, A ,  x);
    else if(alpha == 1)
      MV_Multiply<RangeVector, CrsMatrix, DomainVector, CoeffVector1, CoeffVector2, 1, 2>(betav, y, alphav, A ,  x);
    else if(alpha == -1)
      MV_Multiply < RangeVector, CrsMatrix, DomainVector, CoeffVector1, CoeffVector2, -1, 2 > (betav, y, alphav, A ,  x);
    else
      MV_Multiply<RangeVector, CrsMatrix, DomainVector, CoeffVector1, CoeffVector2, 2, 2>(betav, y, alphav, A ,  x);
  }
}

template < class RangeVector, class CrsMatrix, class DomainVector
, typename Value1, class Layout1, class Device1, class MemoryManagement1
, typename Value2, class Layout2, class Device2, class MemoryManagement2 >
void MV_Multiply(
  const KokkosArray::View<Value1, Layout1, Device1, MemoryManagement1> &betav, const RangeVector &y,
  const KokkosArray::View<Value2, Layout2, Device2, MemoryManagement2> &alphav,
  const CrsMatrix &A , const DomainVector &x)
{
  return MV_Multiply(betav, y, alphav, A, x, 2, 2);
}

template < class RangeVector, class CrsMatrix, class DomainVector
, typename Value1, class Layout1, class Device1, class MemoryManagement1 >
void MV_Multiply(const RangeVector &y,
                 const KokkosArray::View<Value1, Layout1, Device1, MemoryManagement1> &alphav,
                 const CrsMatrix &A , const DomainVector &x)
{
  return MV_Multiply(alphav, y, alphav, A, x, 0, 2);
}

template<class RangeVector, class CrsMatrix, class DomainVector>
void MV_Multiply(const RangeVector &y,
                 const CrsMatrix &A , const DomainVector &x)
{
#ifdef KOKKOS_USE_CUSPARSE
  if(MV_Multiply_Try_CuSparse(0.0, y, 1.0, A, x)) return;
#endif
#ifdef KOKKOS_USE_MKL
  if(MV_Multiply_Try_MKL(0.0, y, 1.0, A, x)) return;
#endif
  typedef KokkosArray::View<typename DomainVector::scalar_type*, typename DomainVector::device_type> aVector;
  aVector a;

  return MV_Multiply(a, y, a, A, x, 0, 1);
}

template<class RangeVector, class CrsMatrix, class DomainVector>
void MV_Multiply(const RangeVector &y, typename DomainVector::const_scalar_type s_a,
                 const CrsMatrix &A , const DomainVector &x)
{
#ifdef KOKKOS_USE_CUSPARSE
  if(MV_Multiply_Try_CuSparse(0.0, y, s_a, A, x)) return;
#endif
#ifdef KOKKOS_USE_MKL
  if(MV_Multiply_Try_MKL(0.0, y, s_a, A, x)) return;
#endif
  typedef KokkosArray::View<typename RangeVector::scalar_type*, typename RangeVector::device_type> aVector;
  aVector a;
  const int numVecs = x.dimension_1();

  if(s_a == -1) {
    return MV_Multiply(a, y, a, A, x, 0, -1);
  } else  if(s_a == 1) {
    return MV_Multiply(a, y, a, A, x, 0, 1);
  }

  if(s_a != 0) {
    a = aVector("a", numVecs);
    typename aVector::HostMirror h_a = KokkosArray::create_mirror_view(a);

    for(int i = 0; i < numVecs; i++)
      h_a(i) = s_a;

    KokkosArray::deep_copy(a, h_a);

    return MV_Multiply(a, y, a, A, x, 0, 2);
  }
}

template<class RangeVector, class CrsMatrix, class DomainVector>
void MV_Multiply(typename RangeVector::const_scalar_type s_b, const RangeVector &y, typename DomainVector::const_scalar_type s_a,
                 const CrsMatrix &A , const DomainVector &x)
{
#ifdef KOKKOS_USE_CUSPARSE
  if(MV_Multiply_Try_CuSparse(s_b, y, s_a, A, x)) return;
#endif
#ifdef KOKKOS_USE_MKL
  if(MV_Multiply_Try_MKL(s_b, y, s_a, A, x)) return;
#endif
  typedef KokkosArray::View<typename RangeVector::scalar_type*, typename RangeVector::device_type> aVector;
  aVector a;
  aVector b;
  int numVecs = x.dimension_1();

    if(s_b == 0) {
      if(s_a == 0)
        return MV_Multiply(a, y, a, A, x, 0, 0);
      else if(s_a == 1)
        return MV_Multiply(a, y, a, A, x, 0, 1);
      else if(s_a == -1)
        return MV_Multiply(a, y, a, A, x, 0, -1);
      else {
        a = aVector("a", numVecs);
        typename aVector::HostMirror h_a = KokkosArray::create_mirror_view(a);

        for(int i = 0; i < numVecs; i++)
          h_a(i) = s_a;

        KokkosArray::deep_copy(a, h_a);
        return MV_Multiply(a, y, a, A, x, 0, 2);
      }
    } else if(s_b == 1) {
      if(s_a == 0)
        return MV_Multiply(a, y, a, A, x, 1, 0);
      else if(s_a == 1)
        return MV_Multiply(a, y, a, A, x, 1, 1);
      else if(s_a == -1)
        return MV_Multiply(a, y, a, A, x, 1, -1);
      else {
        a = aVector("a", numVecs);
        typename aVector::HostMirror h_a = KokkosArray::create_mirror_view(a);

        for(int i = 0; i < numVecs; i++)
          h_a(i) = s_a;

        KokkosArray::deep_copy(a, h_a);
        return MV_Multiply(a, y, a, A, x, 1, 2);
      }
    } else if(s_b == -1) {
      if(s_a == 0)
        return MV_Multiply(a, y, a, A, x, -1, 0);
      else if(s_a == 1)
        return MV_Multiply(a, y, a, A, x, -1, 1);
      else if(s_a == -1)
        return MV_Multiply(a, y, a, A, x, -1, -1);
      else {
        a = aVector("a", numVecs);
        typename aVector::HostMirror h_a = KokkosArray::create_mirror_view(a);

        for(int i = 0; i < numVecs; i++)
          h_a(i) = s_a;

        KokkosArray::deep_copy(a, h_a);
        return MV_Multiply(a, y, a, A, x, -1, 2);
      }
    } else {
      b = aVector("b", numVecs);
      typename aVector::HostMirror h_b = KokkosArray::create_mirror_view(b);

      for(int i = 0; i < numVecs; i++)
        h_b(i) = s_b;

      KokkosArray::deep_copy(b, h_b);

      if(s_a == 0)
        return MV_Multiply(b, y, a, A, x, 2, 0);
      else if(s_a == 1)
        return MV_Multiply(b, y, a, A, x, 2, 1);
      else if(s_a == -1)
        return MV_Multiply(b, y, a, A, x, 2, -1);
      else {
        a = aVector("a", numVecs);
        typename aVector::HostMirror h_a = KokkosArray::create_mirror_view(a);

        for(int i = 0; i < numVecs; i++)
          h_a(i) = s_a;

        KokkosArray::deep_copy(a, h_a);
        return MV_Multiply(b, y, a, A, x, 2, 2);
      }
    }
}

namespace KokkosCrsMatrix
{
template <class CrsMatrixDst, class CrsMatrixSrc>
void deep_copy(CrsMatrixDst A, CrsMatrixSrc B)
{
  KokkosArray::deep_copy(A.graph.entries, B.graph.entries);
  //KokkosArray::deep_copy(A.graph.row_map,B.graph.row_map);
  KokkosArray::deep_copy(A.values, B.values);
}
}

}
#endif /* KOKKOSARRAY_CRSMATRIX_H_ */
