// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __TSQR_Random_NormalGenerator_hpp
#define __TSQR_Random_NormalGenerator_hpp

#include "Tsqr_Impl_Lapack.hpp"
#include <algorithm>
#include <vector>

namespace TSQR {
namespace Random {

/// \class NormalGenerator
/// \brief (Pseudo)random normal(0,1) floating-point number generator.
///
/// \note The implementation currently depends on LAPACK's _LARNV
///   routines.  These are test routines and are not guaranteed to
///   be in the LAPACK library.  They will be if you build LAPACK
///   from source.
template <class Ordinal, class Scalar>
class NormalGenerator {
 private:
  //! Default buffer length.
  static constexpr int defaultBufferLength() { return 100; }

 public:
  typedef Ordinal ordinal_type;
  typedef Scalar scalar_type;

  /// \brief Constructor with custom seed.
  ///
  /// \param iseed [in] Array of four integers, representing the
  ///   seed.  See documentation of LAPACK's _LARNV routines for
  ///   requirements.  In particular, the array elements must be
  ///   in [0,4095], and the last element (iseed[3]) must be odd.
  ///
  /// \param bufferLength [in] How many entries we keep buffered
  ///   at one time.  If you know how many outputs you want, set
  ///   this accordingly, so that all the expense of generation
  ///   happens at construction.
  NormalGenerator(const std::vector<int>& iseed,
                  const int bufferLength = defaultBufferLength())
    : iseed_(4)
    , buffer_(bufferLength)
    , buffer_length_(bufferLength)
    , cur_pos_(0) {
    std::copy(iseed.begin(), iseed.end(), iseed_.begin());
    fill_buffer();
  }

  /// \brief Constructor with default seed.
  ///
  /// The four-integer seed is set to [0, 0, 0, 1], which is a
  /// valid seed and which ensures a reproducible sequence.
  ///
  /// \param bufferLength [in] How many entries we keep buffered at
  ///   one time.  If you know how many outputs you want, set this
  ///   accordingly, so that all the expense of generation happens
  ///   at construction.
  ///
  NormalGenerator(const int bufferLength = defaultBufferLength())
    : iseed_(4)
    , buffer_(bufferLength)
    , buffer_length_(bufferLength)
    , cur_pos_(0) {
    iseed_[0] = 0;
    iseed_[1] = 0;
    iseed_[2] = 0;
    iseed_[3] = 1;
    fill_buffer();
  }

  /// \brief Get the next pseudorandom number.
  ///
  /// If the buffer length is > 0, the buffer is first filled with
  /// newly generated values if it's empty, and then the value is
  /// fetched from the buffer.  Depending on the buffer length,
  /// filling it may take a while, which means that calling this
  /// function may on occasion take longer than you expect.  If
  /// you don't like this behavior, use a buffer length of 1.
  Scalar operator()() { return next(); }

  /// \brief Get the current seed.
  ///
  /// The seed consists of four integers, according to the
  /// requirements of LAPACK's _LARNV routines.  This can be used
  /// to restart the generator, but only if you account for the
  /// buffered values.
  ///
  /// \param iseed [out] Vector of length exactly four.  Resized
  ///   if necessary.
  void
  getSeed(std::vector<int>& iseed) const {
    if (iseed.size() != iseed_.size())
      iseed.resize(iseed_.size());
    std::copy(iseed_.begin(), iseed_.end(), iseed.begin());
  }

 private:
  std::vector<int> iseed_;
  std::vector<Scalar> buffer_;
  int buffer_length_, cur_pos_;

  void
  fill_buffer() {
    Impl::Lapack<Scalar> lapack;

    // LAPACK's _LARNV routine defines this "enum" (just an
    // integer, because it's Fortran) that lets users choose from
    // one of three different pseudorandom distributions:
    // uniform(0,1), uniform(-1,1), and normal(0,1).
    enum distribution_type {
      uniform_0_1  = 1,
      uniform_m1_1 = 2,
      normal_0_1   = 3
    };
    lapack.LARNV(normal_0_1, iseed_.data(), buffer_length_, buffer_.data());
  }

  Scalar
  next() {
    // It's impossible to take the greater-than branch, but we
    // check for robustness' sake.
    if (cur_pos_ >= buffer_length_) {
      fill_buffer();
      cur_pos_ = 0;
    }
    return buffer_[cur_pos_++];
  }
};
}  // namespace Random
}  // namespace TSQR

#endif  // __TSQR_Random_NormalGenerator_hpp
