// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_PACKTRAITS_HPP
#define TPETRA_DETAILS_PACKTRAITS_HPP

///
/// \file Tpetra_Details_PackTraits.hpp
/// \brief Declaration and generic definition of traits class that
///   tells Tpetra::CrsMatrix how to pack and unpack data.
///

#include "Tpetra_ConfigDefs.hpp"
#include "Kokkos_Core.hpp"

namespace Tpetra {
namespace Details {

/// \brief Traits class for packing / unpacking data of type \c T,
///   using Kokkos data structures that live in the given space \c D.
///
/// \tparam T The type of the data to pack / unpack.
/// \tparam D The Kokkos "device" type; where the data live.
template<class T, class D>
struct PackTraits {
  //! The type of data to pack or unpack.
  typedef T value_type;

  /// \brief Whether the number of bytes required to pack one instance
  ///   of \c value_type is fixed at compile time.
  ///
  /// This is true for "plain old data" (POD) types like \c float,
  /// \c double, and \c int.  It is also true of structs or classes of
  /// POD, like Kokkos::complex.  The Sacado and Stokhos packages may
  /// have classes for which this is false.  If false, then the size
  /// of an instance of \c value_type may have been determined at run
  /// time, for example in its constructor.
  ///
  /// Whether or not this is true or false, implementations of
  /// PackTraits may assume that all instances of \c value_type which
  /// the implementation encounters have the same size.
  static const bool compileTimeSize = true;

  //! The type of an input buffer of bytes.
  typedef Kokkos::View<const char*, D, Kokkos::MemoryUnmanaged> input_buffer_type;

  //! The type of an output buffer of bytes.
  typedef Kokkos::View<char*, D, Kokkos::MemoryUnmanaged> output_buffer_type;

  //! The type of an input array of \c value_type.
  typedef Kokkos::View<const value_type*, D, Kokkos::MemoryUnmanaged> input_array_type;

  //! The type of an output array of \c value_type.
  typedef Kokkos::View<value_type*, D, Kokkos::MemoryUnmanaged> output_array_type;

  /// \brief Given an instance of \c value_type allocated with the
  ///   right size, return the "number of values" that make up that
  ///   \c value_type instance.
  ///
  /// This function helps the pack and unpack code that uses
  /// PackTraits correctly handle types that have a size specified at
  /// run time.  PackTraits still assumes that all instances of
  /// \c value_type in an input or output array have the same
  /// run-time size.
  ///
  /// \param x [in] Instance of \c value_type with the correct
  ///   size (possibly determined at run time).
  ///
  /// \return The "number of values" that make up \c x.
  static size_t numValuesPerScalar (const value_type& /* x */) {
    // If your type T is something like Stokhos::UQ::PCE<S>, you must
    // reimplement this function.
    return static_cast<size_t> (1);
  }

  /// \brief Given an instance of \c value_type allocated with the
  ///   right size, allocate and return a one-dimensional array of
  ///   \c value_type.
  ///
  /// This function lets the pack and unpack code that uses PackTraits
  /// correctly handle types that have a size specified at run time.
  /// In particular, it's helpful if that code needs to allocate
  /// temporary buffers of \c value_type.  PackTraits still assumes
  /// that all instances of \c value_type in an input or output array
  /// have the same run-time size.
  ///
  /// \param x [in] Instance of \c value_type with the correct (run-time) size.
  /// \param numEnt [in] Number of entries in the returned array.
  /// \param label [in] Optional string label of the returned
  ///   Kokkos::View.  (Kokkos::View's constructor takes a string
  ///   label, which Kokkos uses for debugging output.)
  ///
  /// \return One-dimensional array of \c value_type, all instances of
  ///   which have the same (run-time) size as \c x.
  ///
  /// \note To implementers of specializations: If the number of bytes
  ///   to pack or unpack your type may be determined at run time, you
  ///   might be able just to use this implementation as-is, and just
  ///   reimplement numValuesPerScalar().
  static Kokkos::View<value_type*, D>
  allocateArray (const value_type& x,
                 const size_t numEnt,
                 const std::string& label = "")
  {
    typedef Kokkos::View<value_type*, D> view_type;
    typedef typename view_type::size_type size_type;

    // This exploits the fact that Kokkos::View's constructor ignores
    // size arguments beyond what the View's type specifies.  For
    // value_type = Stokhos::UQ::PCE<S>, numValuesPerScalar returns
    // something other than 1, and the constructor will actually use
    // that value.
    const size_type numVals = numValuesPerScalar (x);
    return view_type (label, static_cast<size_type> (numEnt), numVals);
  }

  /// \brief Pack the first numEnt entries of the given input buffer
  ///   of \c value_type, into the output buffer of bytes.
  ///
  /// \pre All entries of \c inBuf must have the same (run-time) size.
  ///
  /// \param outBuf [out] Output buffer of bytes (\c char).  Must
  ///   have enough space to hold the packed version of the first
  ///   <tt>numEnt</tt> entries of <tt>inBuf</tt>.
  /// \param inBuf [in] Input buffer of \c value_type.  Must have at
  ///   least \c numEnt entries.
  /// \param numEnt [in] Number of entries to pack.
  ///
  /// \return The number of bytes used to pack \c inBuf into \c outBuf.
  static size_t
  packArray (const output_buffer_type& outBuf,
             const input_array_type& inBuf,
             const size_t numEnt)
  {
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
      // NOTE (mfh 02 Feb 2015) This assumes that all instances of T
      // require the same number of bytes.  To generalize this, we
      // would need to sum up the counts for all entries of inBuf.
      // That of course would suggest that we would need to memcpy
      // each entry separately.
      //
      // We can't just default construct an instance of T, because if
      // T's size is run-time dependent, a default-constructed T might
      // not have the right size.  However, we require that all
      // entries of the input array have the correct size.
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
      memcpy (outBuf.ptr_on_device (), inBuf.ptr_on_device (), numBytes);
      return numBytes;
    }
  }

  /// \brief Unpack \c numEnt \c value_type entries from the given
  ///   input buffer of bytes, to the given output buffer of
  ///   \c value_type.
  ///
  /// \pre All entries of \c outBuf must have the same (run-time)
  ///   size, and that size must be the same as that of the packed
  ///   data that live in \c inBuf.
  ///
  /// \param outBuf [in] Output buffer of \c value_type.  Must have at
  ///   least \c numEnt entries.
  /// \param inBuf [out] Input buffer of bytes (\c char).
  /// \param numEnt [in] Number of \c value_type entries to unpack.
  ///
  /// \return The number of bytes unpacked (i.e., read from \c inBuf).
  static size_t
  unpackArray (const output_array_type& outBuf,
               const input_buffer_type& inBuf,
               const size_t numEnt)
  {
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(
      static_cast<size_t> (outBuf.size ()) < numEnt, std::invalid_argument,
      "PackTraits::unpackArray: outBuf.size() = " << outBuf.size ()
      << " < numEnt = " << numEnt << ".");
#endif // HAVE_TPETRA_DEBUG

    if (numEnt == 0) {
      return static_cast<size_t> (0);
    }
    else {
      // NOTE (mfh 02 Feb 2015) This assumes that all instances of T
      // require the same number of bytes.  To generalize this, we
      // would need to sum up the counts for all entries of inBuf.
      // That of course would suggest that we would need to memcpy
      // each entry separately.
      //
      // We can't just default construct an instance of T, because if
      // T's size is run-time dependent, a default-constructed T might
      // not have the right size.  However, we require that all
      // entries of the input array have the correct size.
      const T& val = packValueCount (outBuf(0));
      const size_t numBytes = numEnt * packValueCount (val);
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION(
        static_cast<size_t> (inBuf.dimension_0 ()) < numBytes,
        std::invalid_argument, "PackTraits::unpackArray: inBuf.dimension_0() = "
        << inBuf.dimension_0 () << " < numBytes = " << numBytes << ".");
#endif // HAVE_TPETRA_DEBUG

      // FIXME (mfh 02,05 Feb 2015) This may assume UVM.  On the other
      // hand, reinterpret_cast may break aliasing and/or alignment
      // rules.
      memcpy (outBuf.ptr_on_device (), inBuf.ptr_on_device (), numBytes);
      return numBytes;
    }
  }

  /// \brief Number of bytes required to pack or unpack the given
  ///   value of type \c value_type.
  ///
  /// \param inVal [in] The input value (see discussion below).
  ///
  /// \return The number of bytes required to pack \c inVal.
  ///
  /// Currently, this function returns the <i>exact</i> amount of
  /// bytes, not an upper bound.  Thus, one can use this function to
  /// predict offsets.  That assumes packing without padding for
  /// (e.g.,) alignment to something bigger than <tt>sizeof(char) =
  /// 1</tt>.  At some point, we may need to extend this to be an
  /// upper bound, rather than an exact value.  Compare to
  /// MPI_PACK_SIZE, which only claims to offer an upper bound.  In
  /// that case, one may not use this function to predict offsets
  /// for unpacking; like MPI_UNPACK, one would need to start at the
  /// beginning of the packed array and unpack sequentially.
  ///
  /// We currently assume that all objects of type \c value_type
  /// require the same number of bytes.  Nevertheless, we require an
  /// instance of \c value_type, in case we want to relax this
  /// assumption in the future.  That's why the brief description of
  /// this function says "the given value of type \c value_type."
  static size_t
  packValueCount (const T& /* inVal */)
  {
    return sizeof (T);
  }

  /// \brief Pack the given value of type \c value_type into the given
  ///   output buffer of bytes (\c char).
  ///
  /// \pre \c outBuf has at least \c packValueCount(inVal) entries.
  ///
  /// \param outBuf [out] Output buffer of bytes.
  /// \param inVal [in] Input value to pack.
  ///
  /// \return The number of bytes used to pack \c inVal.
  static size_t
  packValue (const output_buffer_type& outBuf,
             const T& inVal)
  {
    // It's actually OK for packValueCount to return an upper bound
    // (e.g., padding for alignment).  The memcpy call below will copy
    // any included padding as well as the actual data.
    const size_t numBytes = packValueCount (inVal);

    // FIXME (mfh 02,05 Feb 2015) This may assume UVM.  On the other
    // hand, reinterpret_cast may break aliasing and/or alignment
    // rules.
    memcpy (outBuf.ptr_on_device (), &inVal, numBytes);
    return numBytes;
  }

  /// \brief Unpack the given value from the given output buffer.
  ///
  /// \param outVal [in/out] On output: The unpacked value.
  /// \param inBuf [in] The buffer of packed data from which to unpack
  ///   the output value.
  ///
  /// \return The number of bytes unpacked from \c inBuf.
  ///
  /// We assume that the number of bytes required to pack \c outVal
  /// does not depend on the unpacked data.  That is, \c outVal on
  /// input requires the same number of packed bytes as it should on
  /// output.
  static size_t
  unpackValue (T& outVal, const input_buffer_type& inBuf)
  {
    // It's actually OK for packValueCount to return an upper bound
    // (e.g., padding for alignment).  The memcpy call below will copy
    // any included padding as well as the actual data.
    const size_t numBytes = packValueCount (outVal);

    // FIXME (mfh 02,05 Feb 2015) This may assume UVM.  On the other
    // hand, reinterpret_cast may break aliasing and/or alignment
    // rules.
    memcpy (&outVal, inBuf.ptr_on_device (), numBytes);
    return numBytes;
  }
}; // struct PackTraits

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_PACKTRAITS_HPP
