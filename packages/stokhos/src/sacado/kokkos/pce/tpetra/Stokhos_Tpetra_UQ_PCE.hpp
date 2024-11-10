// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_TPETRA_UQ_PCE_HPP
#define STOKHOS_TPETRA_UQ_PCE_HPP

// This header file should be included whenever compiling any Tpetra
// code with Stokhos scalar types

// MP includes and specializations
#include "Stokhos_Sacado_Kokkos_UQ_PCE.hpp"

// Kokkos includes
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_MultiVector_fwd.hpp"
#include "Tpetra_Vector_fwd.hpp"
#include "Tpetra_Access.hpp"
#include "Kokkos_Core.hpp"
#include <Tpetra_KokkosCompat_ClassicNodeAPI_Wrapper.hpp>
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
      typedef typename std::conditional<
        ::Kokkos::SpaceAccessibility< D, Kokkos::HostSpace>::accessible,
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
      typedef typename std::conditional<
        ::Kokkos::SpaceAccessibility< D, Kokkos::HostSpace>::accessible,
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
#if defined(KOKKOS_ENABLE_SERIAL)
  typedef Kokkos::Serial type;
#else
  typedef Kokkos::HostSpace::execution_space type;
#endif // defined(KOKKOS_ENABLE_SERIAL)
};

template <typename ExecSpace, typename MemSpace>
struct DeviceForNode2< Tpetra::KokkosCompat::KokkosDeviceWrapperNode<ExecSpace, MemSpace> > {
  typedef typename Tpetra::KokkosCompat::KokkosDeviceWrapperNode<ExecSpace, MemSpace>::device_type type;
};

}

#include "Tpetra_Details_PackTraits.hpp"
#include "Tpetra_Details_ScalarViewTraits.hpp"

namespace Tpetra {
namespace Details {

/// \brief Partial specialization of PackTraits for Sacado's PCE UQ type.
///
/// \tparam S The underlying scalar type in the PCE UQ type.
template<class S>
struct PackTraits<Sacado::UQ::PCE<S>> {
  using value_type = Sacado::UQ::PCE<S>;

  /// \brief Whether the number of bytes required to pack one instance
  ///   of \c value_type is fixed at compile time.
  static const bool compileTimeSize = false;

  using input_buffer_type = Kokkos::View<const char*, Kokkos::AnonymousSpace>;
  using output_buffer_type = Kokkos::View<char*, Kokkos::AnonymousSpace>;
  using input_array_type = Kokkos::View<const value_type*, Kokkos::AnonymousSpace>;
  using output_array_type = Kokkos::View<value_type*, Kokkos::AnonymousSpace>;

  using scalar_value_type = typename value_type::value_type;
  using SPT = PackTraits<scalar_value_type>;
  using scalar_input_array_type = typename SPT::input_array_type;
  using scalar_output_array_type = typename SPT::output_array_type;

  KOKKOS_INLINE_FUNCTION
  static size_t numValuesPerScalar (const value_type& x) {
    return x.size ();
  }

  KOKKOS_INLINE_FUNCTION
  static Kokkos::pair<int, size_t>
  packArray (char outBuf[],
             const value_type inBuf[],
             const size_t numEnt)
  {
    using return_type = Kokkos::pair<int, size_t>;
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
    using return_type = Kokkos::pair<int, size_t>;
    size_t numBytes = 0;
    int errorCode = 0;

    if (numEnt == 0) {
      return return_type (errorCode, numBytes);
    }
    else {
      // Check whether output array is contiguously allocated based on the size
      // of the first entry.  We have a simpler method to unpack in this case
      const size_t scalar_size = numValuesPerScalar (outBuf[0]);
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

/// \brief Partial specialization of ScalarViewTraits
///   for Sacado's PCE UQ type.
///
/// \tparam S The underlying scalar type in the PCE UQ type.
/// \tparam D The Kokkos "device" type.
template<typename S, typename D>
struct ScalarViewTraits<Sacado::UQ::PCE<S>, D> {
  using value_type = Sacado::UQ::PCE<S>;
  using device_type = D;

  static Kokkos::View<value_type*, device_type>
  allocateArray (const value_type& x,
                 const size_t numEnt,
                 const std::string& label = "")
  {
    const size_t numVals = PackTraits<value_type>::numValuesPerScalar (x);
    using view_type = Kokkos::View<value_type*, device_type>;
    return view_type (label, numEnt, numVals);
  }
};

} // namespace Details
} // namespace Tpetra

namespace Kokkos {
  template <class S, class L, class G, class N>
  size_t dimension_scalar(const Tpetra::MultiVector<S,L,G,N>& mv) {
    if ( mv.need_sync_device() ) {
      // NOTE (mfh 02 Apr 2019) This doesn't look right.  Shouldn't I
      // want the most recently updated View, which in this case would
      // be the host View?  However, this is what I found when I
      // changed these lines not to call deprecated code, so I'm
      // leaving it.
      return dimension_scalar(mv.getLocalViewHost(Tpetra::Access::ReadOnly));
    }
    return dimension_scalar(mv.getLocalViewDevice(Tpetra::Access::ReadOnly));
  }

  template <class S, class L, class G, class N>
  size_t dimension_scalar(const Tpetra::Vector<S,L,G,N>& v) {
    if ( v.need_sync_device() ) {
      // NOTE (mfh 02 Apr 2019) This doesn't look right.  Shouldn't I
      // want the most recently updated View, which in this case would
      // be the host View?  However, this is what I found when I
      // changed these lines not to call deprecated code, so I'm
      // leaving it.
      return dimension_scalar(v.getLocalViewHost(Tpetra::Access::ReadOnly));
    }
    return dimension_scalar(v.getLocalViewDevice(Tpetra::Access::ReadOnly));
  }
}

#endif // STOKHOS_TPETRA_UQ_PCE_HPP
