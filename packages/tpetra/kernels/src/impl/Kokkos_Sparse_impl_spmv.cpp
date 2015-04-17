#include <Kokkos_Sparse.hpp>

#ifdef KOKKOS_HAVE_OPENMP
#define KOKKOSSPARSE_IMPL_MV_EXEC_SPACE Kokkos::OpenMP
#define KOKKOSSPARSE_IMPL_MV_MEM_SPACE Kokkos::HostSpace
#define KOKKOSSPARSE_IMPL_MV_DEVICE_TYPE Kokkos::Device<KOKKOSSPARSE_IMPL_MV_EXEC_SPACE,KOKKOSSPARSE_IMPL_MV_MEM_SPACE>
#define KOKKOSSPARSE_IMPL_MV_SCALAR double

namespace KokkosSparse {
namespace Impl {

void SPMV<const KOKKOSSPARSE_IMPL_MV_SCALAR, int, KOKKOSSPARSE_IMPL_MV_DEVICE_TYPE, Kokkos::MemoryTraits<Kokkos::Unmanaged>,size_t,
            const KOKKOSSPARSE_IMPL_MV_SCALAR*, Kokkos::LayoutLeft, KOKKOSSPARSE_IMPL_MV_DEVICE_TYPE, Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess>,Kokkos::Impl::ViewDefault,
            KOKKOSSPARSE_IMPL_MV_SCALAR*, Kokkos::LayoutLeft, KOKKOSSPARSE_IMPL_MV_DEVICE_TYPE, Kokkos::MemoryTraits<Kokkos::Unmanaged>,Kokkos::Impl::ViewDefault>::spmv(const char mode[], const Scalar& alpha, const AMatrix& A, const XVector& x, const Scalar& beta, const YVector& y) {
    if( alpha == Kokkos::Details::ArithTraits<Scalar>::zero () ) {
      spmv_alpha<AMatrix,XVector,YVector,0>(mode,alpha,A,x,beta,y);
      return;
    }
    if( alpha == Kokkos::Details::ArithTraits<Scalar>::one () ) {
      spmv_alpha<AMatrix,XVector,YVector,1>(mode,alpha,A,x,beta,y);
      return;
    }
    if( alpha == -Kokkos::Details::ArithTraits<Scalar>::one () ) {
      spmv_alpha<AMatrix,XVector,YVector,-1>(mode,alpha,A,x,beta,y);
      return;
    }
    spmv_alpha<AMatrix,XVector,YVector,2>(mode,alpha,A,x,beta,y);
  }



template<>
struct SPMV<const KOKKOSSPARSE_IMPL_MV_SCALAR, const size_t, KOKKOSSPARSE_IMPL_MV_DEVICE_TYPE, Kokkos::MemoryTraits<Kokkos::Unmanaged>,size_t,
            const KOKKOSSPARSE_IMPL_MV_SCALAR*, Kokkos::LayoutLeft,KOKKOSSPARSE_IMPL_MV_DEVICE_TYPE, Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess>,Kokkos::Impl::ViewDefault,
            KOKKOSSPARSE_IMPL_MV_SCALAR*, Kokkos::LayoutLeft,KOKKOSSPARSE_IMPL_MV_DEVICE_TYPE, Kokkos::MemoryTraits<Kokkos::Unmanaged>,Kokkos::Impl::ViewDefault>
            {
  typedef CrsMatrix<const KOKKOSSPARSE_IMPL_MV_SCALAR, const size_t, KOKKOSSPARSE_IMPL_MV_DEVICE_TYPE, Kokkos::MemoryTraits<Kokkos::Unmanaged>,size_t> AMatrix;
  typedef Kokkos::View<const KOKKOSSPARSE_IMPL_MV_SCALAR*, Kokkos::LayoutLeft,KOKKOSSPARSE_IMPL_MV_DEVICE_TYPE, Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess>,Kokkos::Impl::ViewDefault> XVector;
  typedef Kokkos::View<KOKKOSSPARSE_IMPL_MV_SCALAR*, Kokkos::LayoutLeft,KOKKOSSPARSE_IMPL_MV_DEVICE_TYPE, Kokkos::MemoryTraits<Kokkos::Unmanaged>,Kokkos::Impl::ViewDefault> YVector;
  typedef typename YVector::non_const_value_type Scalar;

  static void spmv(const char mode[], const Scalar& alpha, const AMatrix& A, const XVector& x, const Scalar& beta, const YVector& y) {
    if( alpha == Kokkos::Details::ArithTraits<Scalar>::zero () ) {
      spmv_alpha<AMatrix,XVector,YVector,0>(mode,alpha,A,x,beta,y);
      return;
    }
    if( alpha == Kokkos::Details::ArithTraits<Scalar>::one () ) {
      spmv_alpha<AMatrix,XVector,YVector,1>(mode,alpha,A,x,beta,y);
      return;
    }
    if( alpha == -Kokkos::Details::ArithTraits<Scalar>::one () ) {
      spmv_alpha<AMatrix,XVector,YVector,-1>(mode,alpha,A,x,beta,y);
      return;
    }
    spmv_alpha<AMatrix,XVector,YVector,2>(mode,alpha,A,x,beta,y);
  }

};
}
}
#endif
