// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

/// \file KokkosSparse_SellMatrix.hpp
/// \brief Local sparse sliced ellpack matrix interface
///
/// This file provides KokkosSparse::SellMatrix.  This implements
/// a parse matrix stored in sliced ellpack ("Sell") format.

#ifndef KOKKOSSPARSE_SELLMATRIX_HPP_
#define KOKKOSSPARSE_SELLMATRIX_HPP_

#include <sstream>
#include <stdexcept>

#include "Kokkos_Core.hpp"

namespace KokkosSparse {
namespace Experimental {

/// \class SellMatrix
/// \brief Sliced ell format implementation of a sparse matrix.
/// \tparam ScalarType The type of entries in the sparse matrix.
/// \tparam OrdinalType The type of column indices in the sparse matrix.
/// \tparam Device The Kokkos Device type.
/// \tparam MemoryTraits Traits describing how Kokkos manages and
///   accesses data.  The default parameter suffices for most users.
///
/// "SellMatrix" stands for "sliced ELLPACK" matrix. This format is defined,
/// among other places here: https://doi.org/10.1007/978-3-642-11515-8_10
template <class ScalarType, class OrdinalType, class Device, class MemoryTraits = void,
          class SizeType = KokkosKernels::default_size_type>
class SellMatrix {
  static_assert(std::is_signed<OrdinalType>::value, "SellMatrix requires that OrdinalType is a signed integer type.");

 private:
  using host_mirror_space = typename Kokkos::ViewTraits<ScalarType*, Device, void, MemoryTraits>::host_mirror_space;

 public:
  //! @name Typedefs
  //@{

  //! Type of the matrix's execution space.
  using execution_space = typename Device::execution_space;
  //! Type of the matrix's memory space.
  using memory_space = typename Device::memory_space;
  //! Canonical device type
  using device_type   = Kokkos::Device<execution_space, memory_space>;
  using memory_traits = MemoryTraits;

  //! Type of each value in the matrix.
  using value_type = ScalarType;
  //! Type of each (column) index in the matrix.
  using ordinal_type = OrdinalType;
  //! Type of the offsets used in the track the start and end of each slice.
  using size_type = SizeType;

  //! Type of the view containing the offset to each slice of the matrix.
  using offsets_type = Kokkos::View<size_type*, Kokkos::LayoutRight, device_type, MemoryTraits>;
  //! Const version of the type of the offsets_type in the sparse matrix.
  using const_offset_type = typename offsets_type::const_value_type;
  //! Nonconst version of the type of the entries in the sparse matrix.
  using non_const_offset_type = typename offsets_type::non_const_value_type;
  //! Type of the view containing the column indices of the matrix entries.
  using entries_type = Kokkos::View<ordinal_type*, Kokkos::LayoutRight, device_type, MemoryTraits>;
  //! Const version of the type of the column indices in the sparse matrix.
  using const_entry_type = typename entries_type::const_value_type;
  //! Nonconst version of the type of the column indices in the sparse matrix.
  using non_const_entry_type = typename entries_type::non_const_value_type;
  //! Type of the view containing the values of the matrix entries.
  using values_type = Kokkos::View<value_type*, Kokkos::LayoutRight, device_type, MemoryTraits>;
  //! Const version of the type of the entries in the sparse matrix.
  using const_value_type = typename values_type::const_value_type;
  //! Nonconst version of the type of the entries in the sparse matrix.
  using non_const_value_type = typename values_type::non_const_value_type;

  //! Type of a host-memory mirror of the sparse matrix.
  using host_mirror_type = SellMatrix<ScalarType, OrdinalType, host_mirror_space, MemoryTraits, SizeType>;

  //! Type of a const version of the sparse matrix.
  using const_type = SellMatrix<const_value_type, ordinal_type, device_type, memory_traits, size_type>;

  //@}

  //! Number of rows of the matrix.
  ordinal_type num_rows;
  //! Number of columns of the matrix.
  ordinal_type num_cols;
  //! Number of non-zero values in the matrix.
  size_type nnz;
  //! Number of non-zero in the matrix including padding.
  size_type sell_nnz;
  //! Number of rows per slice.
  ordinal_type num_rows_per_slice;
  //! Number of slices in the matrix.
  ordinal_type num_slices;
  //! Rank-1 view of the offsets for each slice in the sparse matrix.
  offsets_type slice_offsets;
  //! Rank-1 view of column indices of the sparse matrix.
  entries_type entries;
  //! Rank-1 view of values of the sparse matrix.
  values_type values;

  //! @name Constructors and destructor
  //@{

  /// \brief Default constructor that constructs an empty matrix
  KOKKOS_INLINE_FUNCTION SellMatrix() {}

  /// \brief Constructor
  SellMatrix(const OrdinalType nrows, const OrdinalType ncols, const SizeType nnz_, const SizeType padded_nnz,
             const OrdinalType rows_per_slice, const offsets_type& sliceOffsets, const entries_type& cols,
             const values_type& vals)
      : num_rows(nrows),
        num_cols(ncols),
        nnz(nnz_),
        sell_nnz(padded_nnz),
        num_rows_per_slice(rows_per_slice),
        num_slices((nrows + num_rows_per_slice - 1) / num_rows_per_slice),
        slice_offsets(sliceOffsets),
        entries(cols),
        values(vals) {
    // Perform input validation
    if (num_rows_per_slice < num_rows) {  // This requirement should go away once we support proper sell not only ell.
      std::ostringstream os;
      os << "KokkosSparse::SellMatrix: rows_per_slice(" << rows_per_slice << ") < nrows(" << nrows
         << ") which is not supported currently as we only allow one slice that holds the total matrix.";
      throw std::invalid_argument(os.str());
    }
    // Checking that the size of slice_offsets is constistent with the number of slices
    if (num_slices != static_cast<ordinal_type>(slice_offsets.extent(0) - 1)) {
      std::ostringstream os;
      os << "KokkosSparse::SellMatrix: sliceOffsets.extent(0)=" << slice_offsets.extent(0)
         << ", is not equal to num_slices(" << num_slices << ") + 1.";
      throw std::invalid_argument(os.str());
    }
    // The number of values stored in the sell matrix inclues padding so should be larger than nnz
    if (sell_nnz < nnz) {
      std::ostringstream os;
      os << "KokkosSparse::SellMatrix: padded_nz(" << sell_nnz << ") cannot be less than nnz(" << nnz << ").";
      throw std::invalid_argument(os.str());
    }
    // Make sure that extent of the entries and values views match the number of non-zeros in the sell matrix
    if ((static_cast<size_type>(entries.extent(0)) != sell_nnz) ||
        (static_cast<size_type>(values.extent_int(0)) != sell_nnz)) {
      std::ostringstream os;
      os << "KokkosSparse::SellMatrix: entries.extent(0)=" << entries.extent(0)
         << ", values.extent(0)=" << values.extent(0) << " and padded_nnz=" << sell_nnz << " must all be equal.";
      throw std::invalid_argument(os.str());
    }
  }

  //@}
};

/// \class is_sell_matrix
/// \brief is_sell_matrix<T>::value is true if T is a SellMatrix<...>, false
/// otherwise
template <typename>
struct is_sell_matrix : public std::false_type {};
template <typename... P>
struct is_sell_matrix<SellMatrix<P...>> : public std::true_type {};
template <typename... P>
struct is_sell_matrix<const SellMatrix<P...>> : public std::true_type {};

/// \brief Equivalent to is_crs_matrix<T>::value.
template <typename T>
inline constexpr bool is_sell_matrix_v = is_sell_matrix<T>::value;

}  // namespace Experimental
}  // namespace KokkosSparse
#endif  // KOKKOSSPARSE_SELLMATRIX_HPP_
