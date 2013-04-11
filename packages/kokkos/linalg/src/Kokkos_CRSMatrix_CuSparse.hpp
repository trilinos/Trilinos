template<typename T, class RangeVector,class CrsMatrix,class DomainVector>
void MV_Multiply_DoCuSparse(typename KokkosArray::Impl::enable_if<
		!KokkosArray::Impl::is_same<T,double>::value && !KokkosArray::Impl::is_same<T,float>::value, typename RangeVector::scalar_type  >::type s_b
		,const RangeVector & y, typename DomainVector::scalar_type s_a,
		const CrsMatrix & A , const DomainVector & x) {
}

template<typename T, class RangeVector,class CrsMatrix,class DomainVector>
void MV_Multiply_DoCuSparse(typename KokkosArray::Impl::enable_if<KokkosArray::Impl::is_same<T,double>::value, double  >::type s_b
		,const RangeVector & y, double s_a,
		const CrsMatrix & A , const DomainVector & x) {
	cusparseDcsrmm(A.cusparse_handle,CUSPARSE_OPERATION_NON_TRANSPOSE,
	               A.numRows, x.dimension_1(), A.numCols,  A.nnz,
	               &s_a,
	               A.cusparse_descr,
	               A.values.ptr_on_device(),
	               (const int*) A.graph.row_map.ptr_on_device(),
	               A.graph.entries.ptr_on_device(),
	               x.ptr_on_device(),
	               x.dimension_0(),
	               &s_b,
	               y.ptr_on_device(),
				   y.dimension_0());
}

template<typename T, class RangeVector,class CrsMatrix,class DomainVector>
void MV_Multiply_DoCuSparse(typename KokkosArray::Impl::enable_if<KokkosArray::Impl::is_same<T,float>::value, float  >::type s_b
		,const RangeVector & y, float s_a,
		const CrsMatrix & A , const DomainVector & x) {
	cusparseScsrmm(A.cusparse_handle,CUSPARSE_OPERATION_NON_TRANSPOSE,
	               A.numRows, x.dimension_1(), A.numCols,  A.nnz,
	               &s_a,
	               A.cusparse_descr,
	               A.values.ptr_on_device(),
	               (const int*) A.graph.row_map.ptr_on_device(),
	               A.graph.entries.ptr_on_device(),
	               x.ptr_on_device(),
	               x.dimension_0(),
	               &s_b,
	               y.ptr_on_device(),
				   y.dimension_0());
}

//ToDo: strip compatible type attributes (const, volatile); make type of s_b and s_a independent
template<class RangeVector,class CrsMatrix,class DomainVector>
bool MV_Multiply_Try_CuSparse( typename RangeVector::scalar_type s_b,const RangeVector & y, typename DomainVector::scalar_type s_a,
		const CrsMatrix & A , const DomainVector & x)
{
  if(KokkosArray::Impl::is_same<typename RangeVector::scalar_type,float>::value&&
	 KokkosArray::Impl::is_same<typename DomainVector::scalar_type,float>::value&&
	 KokkosArray::Impl::is_same<typename CrsMatrix::scalar_type,float>::value) {
	   MV_Multiply_DoCuSparse<typename RangeVector::scalar_type,RangeVector,CrsMatrix,DomainVector>(s_b,y,s_a,A,x);
	   return true;
  } else
  if(KokkosArray::Impl::is_same<typename RangeVector::scalar_type,double>::value&&
	 KokkosArray::Impl::is_same<typename DomainVector::scalar_type,double>::value&&
	 KokkosArray::Impl::is_same<typename CrsMatrix::scalar_type,double>::value) {
	   MV_Multiply_DoCuSparse<typename RangeVector::scalar_type,RangeVector,CrsMatrix,DomainVector>(s_b,y,s_a,A,x);
	   return true;
  }
  return false;
}
