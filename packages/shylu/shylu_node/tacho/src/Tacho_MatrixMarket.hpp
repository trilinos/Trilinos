// clang-format off
// @HEADER
// *****************************************************************************
//                            Tacho package
//
// Copyright 2022 NTESS and the Tacho contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
// clang-format on
#ifndef __TACHO_MATRIX_MARKET_HPP__
#define __TACHO_MATRIX_MARKET_HPP__

/// \file Tacho_MatrixMarket.hpp
/// \brief IO utilities for matrix market format
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho.hpp"
#include "Tacho_CrsMatrixBase.hpp"

namespace Tacho {

///
/// Coo : Sparse coordinate format; (i, j, val).
///
template <typename ValueType> struct Coo {
  typedef ValueType value_type;

  ordinal_type i, j;
  value_type val;

  Coo() = default;
  Coo(const ordinal_type ii, const ordinal_type jj, const value_type vval) : i(ii), j(jj), val(vval) {}
  Coo(const Coo &b) = default;

  /// \brief Compare "less" index i and j only.
  bool operator<(const Coo &y) const {
    const auto r_val = (this->i - y.i);
    return (r_val == 0 ? this->j < y.j : r_val < 0);
  }

  /// \brief Compare "equality" only index i and j.
  bool operator==(const Coo &y) const { return (this->i == y.i) && (this->j == y.j); }

  /// \brief Compare "in-equality" only index i and j.
  bool operator!=(const Coo &y) const { return !(*this == y); }
};

template <typename T>
inline typename std::enable_if<std::is_same<T, Kokkos::complex<double>>::value ||
                               std::is_same<T, Kokkos::complex<float>>::value>::type
impl_conj_val(T &val, T &conj_val) {
  conj_val = Kokkos::conj(val);
}

template <typename T>
inline typename std::enable_if<std::is_same<T, double>::value || std::is_same<T, float>::value>::type
impl_conj_val(T &val, T &conj_val) {
  conj_val = val;
}

template <typename T>
inline typename std::enable_if<std::is_same<T, Kokkos::complex<double>>::value ||
                               std::is_same<T, Kokkos::complex<float>>::value>::type
impl_is_zero(T &val, bool &is_zero) {
  is_zero = Kokkos::abs(val) < std::numeric_limits<typename T::value_type>::epsilon();
}

template <typename T>
inline typename std::enable_if<std::is_same<T, double>::value || std::is_same<T, float>::value>::type
impl_is_zero(T &val, bool &is_zero) {
  is_zero = std::abs(val) < std::numeric_limits<T>::epsilon();
}

template <typename T>
inline typename std::enable_if<std::is_same<T, double>::value || std::is_same<T, float>::value>::type
impl_read_value_from_file(std::ifstream &file, ordinal_type &row, ordinal_type &col, T &val) {
  file >> row >> col >> val;
}

template <typename T>
inline typename std::enable_if<std::is_same<T, Kokkos::complex<double>>::value ||
                               std::is_same<T, Kokkos::complex<float>>::value>::type
impl_read_value_from_file(std::ifstream &file, ordinal_type &row, ordinal_type &col, T &val) {
  typename T::value_type r, i;
  file >> row >> col >> r >> i;
  val = T(r, i);
}

template <typename T>
inline typename std::enable_if<std::is_same<T, double>::value || std::is_same<T, float>::value>::type
impl_write_value_to_file(std::ofstream &file, const ordinal_type row, const ordinal_type col, const T val,
                         const ordinal_type w = 10) {
  file << std::setw(w) << row << "  " << std::setw(w) << col << "  " << std::setw(w) << val << std::endl;
}

template <typename T>
inline typename std::enable_if<std::is_same<T, Kokkos::complex<double>>::value ||
                               std::is_same<T, Kokkos::complex<float>>::value>::type
impl_write_value_to_file(std::ofstream &file, const ordinal_type row, const ordinal_type col, const T val,
                         const ordinal_type w = 10) {
  file << std::setw(w) << row << "  " << std::setw(w) << col << "  " << std::setw(w) << val.real() << "  "
       << std::setw(w) << val.imag() << std::endl;
}

template <typename ValueType> struct MatrixMarket {

  /// \brief matrix market reader
  template <typename DeviceType>
  static void read(const std::string &filename, CrsMatrixBase<ValueType, DeviceType> &A,
                   const ordinal_type sanitize = 0, const ordinal_type verbose = 0) {
    static_assert(Kokkos::Impl::MemorySpaceAccess<Kokkos::HostSpace, typename DeviceType::memory_space>::assignable,
                  "DeviceType is not assignable from HostSpace");

    Kokkos::Timer timer;

    timer.reset();

    std::ifstream file;
    file.open(filename);

    // reading mm header
    ordinal_type m, n;
    size_type nnz, nnz_input;
    bool symmetry = false, hermitian = false; //, cmplx = false;
    {
      std::string header;
      std::getline(file, header);
      while (file.good()) {
        char c = file.peek();
        if (c == '%' || c == '\n') {
          file.ignore(256, '\n');
          continue;
        }
        break;
      }
      std::transform(header.begin(), header.end(), header.begin(), ::tolower);
      symmetry = (header.find("symmetric") != std::string::npos || header.find("hermitian") != std::string::npos);

      hermitian = (header.find("hermitian") != std::string::npos);

      file >> m >> n >> nnz;
    }

    // read data into coo format
    const ordinal_type mm_base = 1;

    typedef ValueType value_type;
    typedef Coo<value_type> ijv_type;
    std::vector<ijv_type> mm;
    {
      std::vector<ijv_type> mm_org;
      mm_org.reserve(nnz * (symmetry ? 2 : 1));
      for (size_type i = 0; i < nnz; ++i) {
        ordinal_type row, col;
        value_type val;

        impl_read_value_from_file(file, row, col, val);

        row -= mm_base;
        col -= mm_base;

        mm_org.push_back(ijv_type(row, col, val));
        if (symmetry && row != col) {
          value_type conj_val;
          impl_conj_val(val, conj_val);
          mm_org.push_back(ijv_type(col, row, hermitian ? conj_val : val));
        }
      }
      std::sort(mm_org.begin(), mm_org.end(), std::less<ijv_type>());

      // update nnz (this is the nnz from input matrix)
      nnz = mm_org.size();
      nnz_input = nnz;

      // copy to mm
      mm.reserve(nnz);
      if (sanitize) {
        for (size_type i = 0; i < nnz; ++i) {
          bool is_zero;
          impl_is_zero(mm_org[i].val, is_zero);
          if (!is_zero)
            mm.push_back(mm_org[i]);
        }
        nnz = mm.size();
      } else {
        for (size_type i = 0; i < nnz; ++i)
          mm.push_back(mm_org[i]);
        nnz = mm.size();
      }
    }

    // change mm to crs
    Kokkos::View<size_type *, DeviceType> ap("ap", m + 1);
    Kokkos::View<ordinal_type *, DeviceType> aj("aj", nnz);
    Kokkos::View<value_type *, DeviceType> ax("ax", nnz);
    {
      ordinal_type icnt = 0;
      size_type jcnt = 0;
      ijv_type prev = mm[0];

      ap[icnt++] = 0;
      aj[jcnt] = prev.j;
      ax[jcnt++] = prev.val;

      for (auto it = (mm.begin() + 1); it < (mm.end()); ++it) {
        const ijv_type aij = (*it);

        if (aij.i != prev.i)
          ap[icnt++] = jcnt;

        if (aij == prev) {
          aj[jcnt] = aij.j;
          ax[jcnt] += aij.val;
        } else {
          aj[jcnt] = aij.j;
          ax[jcnt++] = aij.val;
        }
        prev = aij;
      }
      ap[icnt++] = jcnt;
      nnz = jcnt;
    }

    // create crs matrix view
    A.clear();
    A.setExternalMatrix(m, n, nnz, ap, aj, ax);

    const double t = timer.seconds();
    if (verbose) {

      printf("Summary: MatrixMarket\n");
      printf("=====================\n");
      printf("  File:      %s\n", filename.c_str());
      printf("  Time\n");
      printf("             time for reading A:                              %10.6f s\n", t);
      printf("\n");
      printf("  Sparse Matrix (%s) \n", (symmetry ? "symmetric" : "non-symmetric"));
      printf("             number of rows:                                  %10d\n", m);
      printf("             number of cols:                                  %10d\n", n);
      printf("             number of nonzeros from input:                   %10d\n", ordinal_type(nnz_input));
      printf("             number of nonzeros after sanitized:              %10d\n", ordinal_type(nnz));
      printf("\n");
    }
  }

  /// \brief matrix marker writer
  template <typename DeviceType>
  static void write(std::ofstream &file, const CrsMatrixBase<ValueType, DeviceType> &A,
                    const int uplo = 0, // 0 - all, 1 - upper, 2 - lower
                    const std::string comment = "%% Tacho::MatrixMarket::Export") {
    static_assert(Kokkos::Impl::MemorySpaceAccess<Kokkos::HostSpace, typename DeviceType::memory_space>::assignable,
                  "DeviceType is not assignable from HostSpace");

    typedef ValueType value_type;
    constexpr bool is_complex = (std::is_same<value_type, Kokkos::complex<float>>::value ||
                                 std::is_same<value_type, Kokkos::complex<double>>::value);

    std::streamsize prec = file.precision();
    file.precision(16);
    file << std::scientific;

    {
      file << "%%MatrixMarket matrix coordinate " << (is_complex ? "complex " : "real ") << std::endl;
      file << comment << std::endl;
    }
    // cnt nnz
    size_type nnz = 0;
    {
      for (ordinal_type i = 0; i < A.NumRows(); ++i) {
        const size_type jbegin = A.RowPtrBegin(i), jend = A.RowPtrEnd(i);
        for (size_type j = jbegin; j < jend; ++j) {
          const auto aj = A.Col(j);
          if (uplo == 1 && i <= aj)
            ++nnz;
          if (uplo == 2 && i >= aj)
            ++nnz;
          if (uplo == 0)
            ++nnz;
        }
      }
      file << A.NumRows() << " " << A.NumCols() << " " << nnz << std::endl;
    }

    const int w = 10;
    {
      for (ordinal_type i = 0; i < A.NumRows(); ++i) {
        const size_type jbegin = A.RowPtrBegin(i), jend = A.RowPtrEnd(i);
        for (size_type j = jbegin; j < jend; ++j) {
          const auto aj = A.Col(j);
          bool flag = false;
          if (uplo == 1 && i <= aj)
            flag = true;
          if (uplo == 2 && i >= aj)
            flag = true;
          if (uplo == 0)
            flag = true;
          if (flag) {
            value_type val = A.Value(j);
            impl_write_value_to_file(file, i + 1, aj + 1, val, w);
          }
        }
      }
    }

    file.unsetf(std::ios::scientific);
    file.precision(prec);
  }
};

} // namespace Tacho

#endif
