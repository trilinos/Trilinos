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
#include "Tpetra_ConfigDefs.hpp"
#include "Kokkos_Core.hpp"
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
      typedef typename Kokkos::Impl::if_c<
        ::Kokkos::Impl::VerifyExecutionCanAccessMemorySpace< D, Kokkos::HostSpace>::value,
        typename D::execution_space, Kokkos::HostSpace>::type
        HostDevice;
      typedef Kokkos::View<T*,D>  view_type;
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
      typedef typename Kokkos::Impl::if_c<
        ::Kokkos::Impl::VerifyExecutionCanAccessMemorySpace< D, Kokkos::HostSpace>::value,
        typename D::execution_space, Kokkos::HostSpace>::type
        HostDevice;
      typedef Kokkos::View<T*,D>  view_type;
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

// Kokkos-Linalg
#include "Kokkos_ArithTraits_UQ_PCE.hpp"
#include "Kokkos_InnerProductSpaceTraits_UQ_PCE.hpp"
#include "Kokkos_MV_UQ_PCE.hpp"
#include "Kokkos_CrsMatrix_UQ_PCE.hpp"
#include "Kokkos_CrsMatrix_UQ_PCE_Cuda.hpp"
#include "Kokkos_TeuchosCommAdapters_UQ_PCE.hpp"
#include "Tpetra_KokkosRefactor_Details_MultiVectorDistObjectKernels_UQ_PCE.hpp"
#include "Tpetra_KokkosRefactor_Details_MultiVectorLocalDeepCopy_UQ_PCE.hpp"
#include "Tpetra_Details_fill_UQ_PCE.hpp"
#include "Kokkos_Random_UQ_PCE.hpp"

namespace Stokhos {

// Traits for determining device type from node type
template <typename Node>
struct DeviceForNode2 {
  // Prefer Serial execution space as the default, but if that's not
  // available, use the Host memory space's default execution space.
#if defined(KOKKOS_HAVE_SERIAL)
  typedef Kokkos::Serial type;
#else
  typedef Kokkos::HostSpace::execution_space type;
#endif // defined(KOKKOS_HAVE_SERIAL)
};

template <typename Device>
struct DeviceForNode2< Kokkos::Compat::KokkosDeviceWrapperNode<Device> > {
  typedef Device type;
};

}

#include "Tpetra_Details_PackTraits.hpp"

namespace Tpetra {
namespace Details {

/// \brief Partial specialization of PackTraits for Sacado's PCE UQ type.
///
/// \tparam S The underlying scalar type in the PCE UQ type.
/// \tparam D The Kokkos "device" type.
template<typename S, typename D>
struct PackTraits< Sacado::UQ::PCE<S>, D > {
  typedef Sacado::UQ::PCE<S> value_type;
  typedef typename D::execution_space execution_space;
  typedef D device_type;
  typedef typename execution_space::size_type size_type;

  /// \brief Whether the number of bytes required to pack one instance
  ///   of \c value_type is fixed at compile time.
  static const bool compileTimeSize = false;

  typedef Kokkos::View<const char*, device_type, Kokkos::MemoryUnmanaged> input_buffer_type;
  typedef Kokkos::View<char*, device_type, Kokkos::MemoryUnmanaged> output_buffer_type;
  typedef Kokkos::View<const value_type*, device_type, Kokkos::MemoryUnmanaged> input_array_type;
  typedef Kokkos::View<value_type*, device_type, Kokkos::MemoryUnmanaged> output_array_type;

  typedef typename value_type::value_type scalar_value_type;
  typedef PackTraits< scalar_value_type, device_type > SPT;
  typedef typename SPT::input_array_type scalar_input_array_type;
  typedef typename SPT::output_array_type scalar_output_array_type;

  KOKKOS_INLINE_FUNCTION
  static size_t numValuesPerScalar (const value_type& x) {
    return x.size ();
  }

  static Kokkos::View<value_type*, device_type>
  allocateArray (const value_type& x, const size_t numEnt, const std::string& label = "")
  {
    typedef Kokkos::View<value_type*, device_type> view_type;

    const size_type numVals = numValuesPerScalar (x);
    return view_type (label, static_cast<size_type> (numEnt), numVals);
  }

  KOKKOS_INLINE_FUNCTION
  static Kokkos::pair<int, size_t>
  packArray (char outBuf[],
             const value_type inBuf[],
             const size_t numEnt)
  {
    typedef Kokkos::pair<int, size_t> return_type;
    size_t numBytes = 0;
    int errorCode = 0;

    if (numEnt == 0) {
      return return_type (errorCode, numBytes);
    }
    else {
      // Check whether input array is contiguously allocated based on the size
      // of the first entry.  We can only pack contiguously allocated data
      // since that is the only way we can guarrantee all of the PCE arrays
      // are the same size and the buffer will allocated correctly.
      const size_t scalar_size = numValuesPerScalar (inBuf[0]);
      const scalar_value_type* last_coeff = inBuf[numEnt - 1].coeff ();
      const scalar_value_type* last_coeff_expected =
        inBuf[0].coeff () + (numEnt - 1)*scalar_size;
      const bool is_contiguous = (last_coeff == last_coeff_expected);
      if (! is_contiguous) {
        // Cannot pack non-contiguous PCE array since buffer size calculation
        // is likely wrong.
        errorCode = 3;
        return return_type (errorCode, numBytes);
      }

      // Check we are packing length-1 PCE arrays (mean-based preconditioner).
      // We can technically pack length > 1, but the unpack assumes the
      // output array is sized appropriately.  Currently this is not the case
      // in Tpetra::CrsMatrix::transferAndFillComplete() which allocates a
      // local Teuchos::Array for the CSR values, which will only be length-1
      // by default.
      if (scalar_size != 1) {
        // Cannot pack PCE array with pce_size > 1 since unpack array
        // may not be allocated correctly.
        errorCode = 4;
        return return_type (errorCode, numBytes);
      }

      const size_t flat_numEnt = numEnt * scalar_size;
      return SPT::packArray (outBuf, inBuf[0].coeff (), flat_numEnt);
    }
  }

  KOKKOS_INLINE_FUNCTION
  static Kokkos::pair<int, size_t>
  unpackArray (value_type outBuf[],
               const char inBuf[],
               const size_t numEnt)
  {
    typedef Kokkos::pair<int, size_t> return_type;
    size_t numBytes = 0;
    int errorCode = 0;

    if (numEnt == 0) {
      return return_type (errorCode, numBytes);
    }
    else {
      // Check whether output array is contiguously allocated based on the size
      // of the first entry.  We have a simpler method to unpack in this case
      const size_type scalar_size = numValuesPerScalar (outBuf[0]);
      const scalar_value_type* last_coeff = outBuf[numEnt - 1].coeff ();
      const scalar_value_type* last_coeff_expected =
        outBuf[0].coeff () + (numEnt - 1) * scalar_size;
      const bool is_contiguous = (last_coeff == last_coeff_expected);

      if (is_contiguous) {
        // Unpack all of the PCE coefficients for the whole array
        const size_t flat_numEnt = numEnt * scalar_size;
        return SPT::unpackArray (outBuf[0].coeff (), inBuf, flat_numEnt);
      }
      else {
        // Unpack one entry at a time.  This assumes each entry of outBuf
        // is the correct size based on the packing.  This is is only
        // guarranteed to be true for pce_size == 1, hence the check in
        // packArray().
        size_t numBytesTotal = 0;
        for (size_t i = 0; i < numEnt; ++i) {
          const char* inBufVal = inBuf + numBytesTotal;
          const size_t numBytes = unpackValue (outBuf[i], inBufVal);
          numBytesTotal += numBytes;
        }
        return return_type (errorCode, numBytesTotal);
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  static size_t
  packValueCount (const value_type& inVal)
  {
    return inVal.size () * SPT::packValueCount (inVal.val ());
  }

  KOKKOS_INLINE_FUNCTION
  static size_t
  packValue (char outBuf[],
             const value_type& inVal)
  {
    const size_t numBytes = packValueCount (inVal);
    memcpy (outBuf, inVal.coeff (), numBytes);
    return numBytes;
  }

  KOKKOS_INLINE_FUNCTION
  static size_t
  packValue (char outBuf[],
             const size_t outBufIndex,
             const value_type& inVal)
  {
    const size_t numBytes = packValueCount (inVal);
    const size_t offset = outBufIndex * numBytes;
    memcpy (outBuf + offset, inVal.coeff (), numBytes);
    return numBytes;
  }

  KOKKOS_INLINE_FUNCTION
  static size_t
  unpackValue (value_type& outVal, const char inBuf[])
  {
    const size_t numBytes = packValueCount (outVal);
    memcpy (outVal.coeff (), inBuf, numBytes);
    return numBytes;
  }
}; // struct PackTraits

} // namespace Details
} // namespace Tpetra

namespace Tpetra {
  template <class S, class L, class G, class N> class MultiVector;
  template <class S, class L, class G, class N> class Vector;
}

namespace Kokkos {
  template <class S, class L, class G, class N>
  size_t dimension_scalar(const Tpetra::MultiVector<S,L,G,N>& mv) {
    typedef Tpetra::MultiVector<S,L,G,N> MV;
    typedef typename MV::dual_view_type dual_view_type;
    typedef typename dual_view_type::t_dev device_type;
    typedef typename dual_view_type::t_host host_type;
    dual_view_type dual_view = mv.getDualView();
    if (dual_view.modified_host() > dual_view.modified_device())
      return dimension_scalar(dual_view.template view<device_type>());
    return dimension_scalar(dual_view.template view<host_type>());
  }

  template <class S, class L, class G, class N>
  size_t dimension_scalar(const Tpetra::Vector<S,L,G,N>& v) {
    typedef Tpetra::Vector<S,L,G,N> V;
    typedef typename V::dual_view_type dual_view_type;
    typedef typename dual_view_type::t_dev device_type;
    typedef typename dual_view_type::t_host host_type;
    dual_view_type dual_view = v.getDualView();
    if (dual_view.modified_host() > dual_view.modified_device())
      return dimension_scalar(dual_view.template view<device_type>());
    return dimension_scalar(dual_view.template view<host_type>());
  }
}

#endif // STOKHOS_TPETRA_UQ_PCE_HPP
