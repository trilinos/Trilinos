// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

/// \brief Traits class for packing / unpacking data of type \c T.
///
/// \tparam T The type of the data to pack / unpack.
template<class T>
struct PackTraits {
  //! The type of data to pack or unpack.
  using value_type = T;

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
  using input_buffer_type = Kokkos::View<const char*, Kokkos::AnonymousSpace>;

  //! The type of an output buffer of bytes.
  using output_buffer_type = Kokkos::View<char*, Kokkos::AnonymousSpace>;

  //! The type of an input array of \c value_type.
  using input_array_type = Kokkos::View<const value_type*, Kokkos::AnonymousSpace>;

  //! The type of an output array of \c value_type.
  using output_array_type = Kokkos::View<value_type*, Kokkos::AnonymousSpace>;

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
  KOKKOS_INLINE_FUNCTION
  static size_t
  numValuesPerScalar (const value_type& /* x */) {
    // If your type T is something like Stokhos::UQ::PCE<S>, you must
    // reimplement this function.
    return static_cast<size_t> (1);
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
  /// \return {The number of bytes used to pack \c inBuf into \c outBuf,
  //    packArray error code (success = 0)}
  KOKKOS_INLINE_FUNCTION
  static Kokkos::pair<int, size_t>
  packArray (char outBuf[],
             const value_type inBuf[],
             const size_t numEnt)
  {
    size_t numBytes = 0;
    int errorCode = 0;
    typedef Kokkos::pair<int, size_t> pair_type;

    if (numEnt == 0) {
      return pair_type (errorCode, numBytes);
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
      numBytes = numEnt * packValueCount (inBuf[0]);

      // As of CUDA 6, it's totally fine to use memcpy in a CUDA device
      // function.  It does what one would expect.
      memcpy (outBuf, inBuf, numBytes);
      return pair_type (errorCode, numBytes);
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
  /// \return {The number of bytes unpacked (i.e., read from \c inBuf).
  //    unpackArray error code (success = 0)}
  KOKKOS_INLINE_FUNCTION
  static Kokkos::pair<int, size_t>
  unpackArray (value_type outBuf[],
               const char inBuf[],
               const size_t numEnt)
  {
    size_t numBytes = 0;
    int errorCode = 0;
    typedef Kokkos::pair<int, size_t> pair_type;

    if (numEnt == 0) {
      return pair_type (errorCode, numBytes);
    }
    else {
      // NOTE (mfh 02 Feb 2015) This assumes that all instances of
      // value_type require the same number of bytes.  To generalize
      // this, we would need to sum up the counts for all entries of
      // inBuf.  That of course would suggest that we would need to
      // memcpy each entry separately.
      //
      // We can't just default construct an instance of value_type,
      // because if value_type's size is run-time dependent, a
      // default-constructed value_type might not have the right size.
      // However, we require that all entries of the input array have
      // the correct size.
      const value_type& val = outBuf[0];
      numBytes = numEnt * packValueCount (val);

      // As of CUDA 6, it's totally fine to use memcpy in a CUDA device
      // function.  It does what one would expect.
      memcpy ((void*) outBuf, inBuf, numBytes);
      return pair_type (errorCode, numBytes);
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
  KOKKOS_INLINE_FUNCTION
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
  KOKKOS_INLINE_FUNCTION
  static size_t
  packValue (char outBuf[],
             const T& inVal)
  {
    // It's actually OK for packValueCount to return an upper bound
    // (e.g., padding for alignment).  The memcpy call below will copy
    // any included padding as well as the actual data.
    const size_t numBytes = packValueCount (inVal);

    // As of CUDA 6, it's totally fine to use memcpy in a CUDA device
    // function.  It does what one would expect.
    memcpy (outBuf, &inVal, numBytes);
    return numBytes;
  }

  /// \brief Pack the given value of type \c value_type into the given
  ///   output buffer of bytes (\c char).
  ///
  /// \pre \c outBuf has at least \c packValueCount(inVal) entries.
  ///
  /// \param outBuf [out] Output buffer of bytes.
  /// \param outBufIndex [in] Index into output buffer (multiplied by
  ///   the number of bytes needed to pack inVal).
  /// \param inVal [in] Input value to pack.
  ///
  /// \return The number of bytes used to pack \c inVal.
  KOKKOS_INLINE_FUNCTION
  static size_t
  packValue (char outBuf[],
             const size_t outBufIndex,
             const T& inVal)
  {
    // It's actually OK for packValueCount to return an upper bound
    // (e.g., padding for alignment).  The memcpy call below will copy
    // any included padding as well as the actual data.
    const size_t numBytes = packValueCount (inVal);
    const size_t offset = outBufIndex * numBytes;

    // As of CUDA 6, it's totally fine to use memcpy in a CUDA device
    // function.  It does what one would expect.
    memcpy (outBuf + offset, &inVal, numBytes);
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
  KOKKOS_INLINE_FUNCTION
  static size_t
  unpackValue (T& outVal, const char inBuf[])
  {
    // It's actually OK for packValueCount to return an upper bound
    // (e.g., padding for alignment).  The memcpy call below will copy
    // any included padding as well as the actual data.
    const size_t numBytes = packValueCount (outVal);

    // As of CUDA 6, it's totally fine to use memcpy in a CUDA device
    // function.  It does what one would expect.
    memcpy ((void*) &outVal, inBuf, numBytes);
    return numBytes;
  }
}; // struct PackTraits

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_PACKTRAITS_HPP
