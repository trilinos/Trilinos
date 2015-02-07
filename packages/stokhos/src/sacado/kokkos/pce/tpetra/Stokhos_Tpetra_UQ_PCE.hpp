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

// Kokkos-Linalg
#include "Kokkos_ArithTraits_UQ_PCE.hpp"
#include "Kokkos_InnerProductSpaceTraits_UQ_PCE.hpp"
#include "Kokkos_MV_UQ_PCE.hpp"
#include "Kokkos_CrsMatrix_UQ_PCE.hpp"
#include "Kokkos_CrsMatrix_UQ_PCE_Cuda.hpp"
#include "Kokkos_TeuchosCommAdapters_UQ_PCE.hpp"
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

#include "Tpetra_Import_Util2.hpp"

namespace Tpetra {
namespace Import_Util {

/// \brief Partial specialization of PackTraits for Sacado's PCE UQ type.
///
/// \tparam S The underlying scalar type in the PCE UQ type.
/// \tparam D The Kokkos "device" type.
template<typename S, typename D>
struct PackTraits< Sacado::UQ::PCE<S>, D > {
  typedef Sacado::UQ::PCE<S> value_type;

  /// \brief Whether the number of bytes required to pack one instance
  ///   of \c value_type is fixed at compile time.
  static const bool compileTimeSize = false;

  typedef Kokkos::View<const char*, D, Kokkos::MemoryUnmanaged> input_buffer_type;
  typedef Kokkos::View<char*, D, Kokkos::MemoryUnmanaged> output_buffer_type;
  typedef Kokkos::View<const value_type*, D, Kokkos::MemoryUnmanaged> input_array_type;
  typedef Kokkos::View<value_type*, D, Kokkos::MemoryUnmanaged> output_array_type;

  static size_t numValuesPerScalar (const S& x) {
    return x.size ();
  }

  static Kokkos::View<S*, D>
  allocateArray (const S& x, const size_t numEnt, const std::string& label = "")
  {
    typedef Kokkos::View<S*, D> view_type;
    typedef typename view_type::size_type size_type;

    const size_type numVals = numValuesPerScalar (x);
    return view_type (label, static_cast<size_type> (numEnt), numVals);
  }

  static size_t
  packArray (const output_buffer_type& outBuf,
             const input_array_type& inBuf,
             const size_t numEnt)
  {
    typedef typename value_type::value_type SVT; // "scalar value type"

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(
      static_cast<size_t> (inBuf.dimension_0 ()) < numEnt,
      std::invalid_argument, "PackTraits::packArray: inBuf.dimension_0() = "
      << inBuf.dimension_0 () << " < numEnt = " << numEnt << ".");
#endif // HAVE_TPETRA_DEBUG

    if (numEnt == 0) {
      return 0;
    }
    else {
      // NOTE (mfh 02 Feb 2015) This assumes that all instances of
      // value_type require the same number of bytes.  To generalize
      // this, we would need to sum up the counts for all entries of
      // inBuf.  That of course would suggest that we would need to
      // memcpy each entry separately.
      //
      // We can't just default construct an instance of value_type,
      // because value_type's size is run-time dependent.  However, we
      // require that all entries of the input array have the correct
      // size, so it suffices to ask the first entry of the input
      // array for its size.
      const size_t numBytes = numEnt * packValueCount (inBuf(0));
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION(
        static_cast<size_t> (outBuf.dimension_0 ()) < numBytes,
        std::invalid_argument, "PackTraits::packArray: outBuf.dimension_0() = "
        << outBuf.dimension_0 () << " < numBytes = " << numBytes << ".");
#endif // HAVE_TPETRA_DEBUG

      // FIXME (mfh 02,05 Feb 2015) This may assume UVM.  On the other
      // hand, reinterpret_cast may break aliasing and/or alignment
      // rules.
      const SVT* inBufRaw = inBuf(0).coeff ();
      memcpy (outBuf.ptr_on_device (), inBufRaw, numBytes);
      return numBytes;
    }
  }

  static size_t
  unpackArray (const output_array_type& outBuf,
               const input_buffer_type& inBuf,
               const size_t numEnt)
  {
    typedef typename value_type::value_type SVT; // "scalar value type"

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(
      static_cast<size_t> (outBuf.dimension_0 ()) < numEnt, std::invalid_argument,
      "PackTraits::unpackArray: outBuf.dimension_0 () = " << outBuf.dimension_0 ()
      << " < numEnt = " << numEnt << ".");
#endif // HAVE_TPETRA_DEBUG

    if (numEnt == 0) {
      return static_cast<size_t> (0);
    }
    else {
      // NOTE (mfh 02 Feb 2015) This assumes that all instances of
      // value_type require the same number of bytes.  To generalize
      // this, we would need to sum up the counts for all entries of
      // outBuf.  That of course would suggest that we would need to
      // memcpy each entry separately.
      //
      // We can't just default construct an instance of value_type,
      // because if value_type's size is run-time dependent, a
      // default-constructed value_type might not have the right size.
      // However, we require that all entries of the output array have
      // the correct size, so it suffices to look at the first entry.
      const size_t numBytes = numEnt * packValueCount (outBuf(0));
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION(
        static_cast<size_t> (inBuf.dimension_0 ()) < numBytes,
        std::invalid_argument, "PackTraits::unpackArray: inBuf.dimension_0() = "
        << inBuf.dimension_0 () << " < numBytes = " << numBytes << ".");
#endif // HAVE_TPETRA_DEBUG

      // FIXME (mfh 02,05 Feb 2015) This may assume UVM.  On the other
      // hand, reinterpret_cast may break aliasing and/or alignment
      // rules.
      SVT* outBufRaw = outBuf(0).coeff ();
      memcpy (outBufRaw, inBuf.ptr_on_device (), numBytes);
      return numBytes;
    }
  }

  static size_t
  packValueCount (const value_type& inVal)
  {
    typedef typename value_type::value_type SVT; // "scalar value type"
    // NOTE (mfh 06 Feb 2015) It might be reasonable to assume that
    // SVT is default constructible, and that all SVT instances have
    // the same size.  On the other hand, the latter might not be true
    // if Stokhos allows nesting of PCE types.  It's safer just to ask
    // inVal for its (zeroth) value.
    return inVal.size () * PackTraits<SVT, D>::packValueCount (inVal.val ());
  }

  static size_t
  packValue (const output_buffer_type& outBuf,
             const value_type& inVal)
  {
    const size_t numBytes = packValueCount (inVal);
    memcpy (outBuf.ptr_on_device (), inVal.coeff (), numBytes);
    return numBytes;
  }

  static size_t
  unpackValue (value_type& outVal, const input_buffer_type& inBuf)
  {
    const size_t numBytes = packValueCount (outVal);
    memcpy (outVal.coeff (), inBuf.ptr_on_device (), numBytes);
    return numBytes;
  }
}; // struct PackTraits

} // namespace Import_Util
} // namespace Tpetra

#endif // STOKHOS_TPETRA_UQ_PCE_HPP
