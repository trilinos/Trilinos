// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_TPETRA_UQ_PCE_HPP
#define STOKHOS_TPETRA_UQ_PCE_HPP

// This header file should be included whenever compiling any Tpetra
// code with Stokhos scalar types

// MP includes and specializations
#include "Stokhos_Sacado_Kokkos_UQ_PCE.hpp"

// Kokkos includes
#include "KokkosClassic_config.h"
#include "Kokkos_Core.hpp"
#if defined(HAVE_KOKKOSCLASSIC_KOKKOSCOMPAT)
#include "Kokkos_BufferMacros.hpp"
#include "KokkosCompat_ClassicNodeAPI_Wrapper.hpp"
#include "KokkosCompat_View.hpp"
#include "KokkosCompat_View_def.hpp"

// Hack for mean-based prec-setup where we get the PCE size from the first
// entry in the Teuchos::ArrayView.  This is clearly quite dangerous and is
// likely to bite us in the ass at some point!
namespace Kokkos {
  namespace Compat {
    template <typename D, typename S>
    Kokkos::View<Sacado::UQ::PCE<S>*,D>
    getKokkosViewDeepCopy(const Teuchos::ArrayView< Sacado::UQ::PCE<S> >& a) {
      typedef Sacado::UQ::PCE<S> T;
      typedef Kokkos::View<T*,D>  view_type;
      typedef typename view_type::host_mirror_space HostDevice;
      typedef Kokkos::View<T*,typename view_type::array_layout,HostDevice,Kokkos::MemoryUnmanaged> unmanaged_host_view_type;
      if (a.size() == 0)
        return view_type();
      view_type v("", a.size(), a[0].size());
      unmanaged_host_view_type hv(a.getRawPtr(), a.size(), a[0].size());
      Kokkos::deep_copy(v,hv);
      return v;
    }

    template <typename D, typename S>
    Kokkos::View<const Sacado::UQ::PCE<S>*,D>
    getKokkosViewDeepCopy(const Teuchos::ArrayView<const Sacado::UQ::PCE<S> >& a) {
      typedef Sacado::UQ::PCE<S> T;
      typedef Kokkos::View<T*,D>  view_type;
      typedef typename view_type::host_mirror_space HostDevice;
      typedef Kokkos::View<const T*,typename view_type::array_layout,HostDevice,Kokkos::MemoryUnmanaged> unmanaged_host_view_type;
      if (a.size() == 0)
        return view_type();
      view_type v("", a.size(), a[0].size());
      unmanaged_host_view_type hv(a.getRawPtr(), a.size(), a[0].size());
      Kokkos::deep_copy(v,hv);
      return v;
    }
  }
}
#endif

// Kokkos-Linalg
#include "Tpetra_ConfigDefs.hpp"
#if defined(TPETRA_HAVE_KOKKOS_REFACTOR)
#include "Kokkos_ArithTraits_UQ_PCE.hpp"
#include "Kokkos_InnerProductSpaceTraits_UQ_PCE.hpp"
#include "Kokkos_MV_UQ_PCE.hpp"
#include "Kokkos_CrsMatrix_UQ_PCE.hpp"
#include "Kokkos_CrsMatrix_UQ_PCE_Cuda.hpp"
#include "Kokkos_TeuchosCommAdapters_UQ_PCE.hpp"
#include "Kokkos_Random_UQ_PCE.hpp"
#endif

namespace Stokhos {

// Traits for determining device type from node type
template <typename Node>
struct DeviceForNode2 {
  typedef Kokkos::Serial type;
};

#if defined(HAVE_KOKKOSCLASSIC_KOKKOSCOMPAT)
template <typename Device>
struct DeviceForNode2< Kokkos::Compat::KokkosDeviceWrapperNode<Device> > {
  typedef Device type;
};
#endif

}

#if defined(TPETRA_HAVE_KOKKOS_REFACTOR)
#include "Tpetra_Import_Util2.hpp"
namespace Tpetra {
  namespace Import_Util {
    template <typename S, typename LO, typename GO, typename D>
    struct MatrixSerializationTraits<
      CrsMatrix< Sacado::UQ::PCE<S>,LO,GO,Kokkos::Compat::KokkosDeviceWrapperNode<D> > > {
      typedef Sacado::UQ::PCE<S> Scalar;
      typedef Kokkos::Compat::KokkosDeviceWrapperNode<D> Node;
      typedef CrsMatrix<Scalar,LO,GO,Node> Matrix;

      typedef typename Scalar::value_type scalar_value;

      static inline
      size_t scalarSize( const Matrix& mat ) {
        const size_t pce_size = mat.getLocalMatrix().values.sacado_size();
        return pce_size *sizeof(scalar_value);
      }

      static inline
      void packBuffer( const Matrix& mat,
                       const size_t numEntries,
                       const Teuchos::ArrayView<const Scalar>& vals,
                       const Teuchos::ArrayView<char> packed_vals ) {
        if (numEntries == 0) return;
        const size_t pce_size = mat.getLocalMatrix().values.sacado_size();
        const scalar_value* pce_vals_beg = vals[0].coeff();
        const scalar_value* pce_vals_end = vals[numEntries-1].coeff()+pce_size;
        const bool is_contiguous =
          ( pce_vals_end == pce_vals_beg + numEntries*pce_size );
        TEUCHOS_TEST_FOR_EXCEPTION( !is_contiguous, std::logic_error,
                                    "PCE array is not contiguous!" );
        scalar_value* packed_vals_scalar =
          reinterpret_cast<scalar_value*>(packed_vals.getRawPtr());
        std::copy( pce_vals_beg, pce_vals_end, packed_vals_scalar );
      }

      static inline
      void unpackScalar( const Matrix& mat,
                         const char * val_char,
                         Scalar& val ) {
        const size_t pce_size = mat.getLocalMatrix().values.sacado_size();
        scalar_value* pce_vals =
          const_cast<scalar_value*>(reinterpret_cast<const scalar_value*>(val_char));
        val = Scalar( mat.getLocalMatrix().values.cijk(), pce_size, pce_vals,
                      false );
      }
    };
  }
}

#endif

#endif // STOKHOS_TPETRA_UQ_PCE_HPP
