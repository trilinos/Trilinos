//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <vector>
#include <unordered_set>
#include <stdexcept>
#include <type_traits>
#ifndef _KOKKOSKERNELSIOUTILS_HPP
#define _KOKKOSKERNELSIOUTILS_HPP

#include "Kokkos_ArithTraits.hpp"
#include <Kokkos_Core.hpp>
#include "Kokkos_Random.hpp"
#include "KokkosKernels_SimpleUtils.hpp"
#include <sys/stat.h>

namespace KokkosKernels {

namespace Impl {

// Get the interval for Kokkos::fill_random
// For real, interval is (-mag, mag)
// For complex, both real and imaginary parts will have interval (-mag, mag)
template <typename Scalar>
inline void getRandomBounds(double mag, Scalar &start, Scalar &end) {
  start = -mag * Kokkos::ArithTraits<Scalar>::one();
  end   = mag * Kokkos::ArithTraits<Scalar>::one();
}

template <>
inline void getRandomBounds(double mag, Kokkos::complex<float> &start, Kokkos::complex<float> &end) {
  start = Kokkos::complex<float>(-mag, -mag);
  end   = Kokkos::complex<float>(mag, mag);
}

template <>
inline void getRandomBounds(double mag, Kokkos::complex<double> &start, Kokkos::complex<double> &end) {
  start = Kokkos::complex<double>(-mag, -mag);
  end   = Kokkos::complex<double>(mag, mag);
}

template <typename stype>
void md_malloc(stype **arr, size_t n, std::string /*alloc_str*/ = "") {
  *arr = new stype[n];
  if (*arr == NULL) {
    throw std::runtime_error("Memory Allocation Problem\n");
  }
}

template <typename idx, typename wt>
struct Edge {
  idx src;
  idx dst;
  wt ew;
  bool operator<(const Edge<idx, wt> &a) const {
    // return !((this->src < a.src) || (this->src == a.src && this->dst <
    // a.dst));
    return (this->src < a.src) || (this->src == a.src && this->dst < a.dst);
  }
};

////////////////////////////////////////////////////////////////////////////////
// From MTGL
////////////////////////////////////////////////////////////////////////////////
inline size_t kk_get_file_size(const char *file) {
  // struct stat stat_buf;

#ifdef _WIN32
  struct _stat stat_buf;
  int retval = _stat(file, &stat_buf);
#else
  struct stat stat_buf;
  int retval = stat(file, &stat_buf);
#endif

  return retval == 0 ? size_t(stat_buf.st_size) : size_t(0);
}

template <typename lno_t>
void buildEdgeListFromBinSrcTarg_undirected(const char *fnameSrc, const char *fnameTarg, size_t &numEdges, lno_t **srcs,
                                            lno_t **dst) {
  size_t srcFileSize = kk_get_file_size(fnameSrc);
  size_t trgFileSize = kk_get_file_size(fnameTarg);
  // test these values

  size_t srcSize = srcFileSize / sizeof(lno_t);
  size_t trgSize = trgFileSize / sizeof(lno_t);
  if (srcSize != trgSize) {
    throw std::runtime_error("Src and Target file needs to be the same size");
  }
  // Assumption that each edge is listed once
  numEdges = srcSize;

  md_malloc<lno_t>(srcs, numEdges);
  md_malloc<lno_t>(dst, numEdges);

  ////////////////////////////////////////////////////////
  // Read source data into buffer
  ////////////////////////////////////////////////////////

  std::ifstream myFile(fnameSrc, std::ios::in | std::ios::binary);

  myFile.read((char *)*srcs, sizeof(lno_t) * (numEdges));

  myFile.close();

  std::ifstream myFile2(fnameTarg, std::ios::in | std::ios::binary);

  myFile2.read((char *)*dst, sizeof(lno_t) * (numEdges));

  myFile2.close();
  //
}

template <typename idx_array_type>
inline void kk_write_1Dview_to_file(idx_array_type view, const char *filename) {
  typedef typename idx_array_type::HostMirror host_type;
  // typedef typename idx_array_type::size_type idx;
  host_type host_view = Kokkos::create_mirror_view(view);
  Kokkos::deep_copy(host_view, view);
  Kokkos::fence();
  std::ofstream myFile(filename, std::ios::out);
  for (size_t i = 0; i < view.extent(0); ++i) {
    myFile << host_view(i) << std::endl;
  }
  myFile.close();
}

template <typename idx_array_type>
inline void kk_read_1Dview_from_file(idx_array_type &view, const char *filename) {
  typedef typename idx_array_type::HostMirror host_type;
  // typedef typename idx_array_type::size_type idx;
  host_type host_view = Kokkos::create_mirror_view(view);
  std::ifstream myFile(filename, std::ios::in);

  for (size_t i = 0; i < view.extent(0); ++i) {
    myFile >> host_view(i);
  }
  myFile.close();
  Kokkos::deep_copy(view, host_view);
  Kokkos::fence();
}

template <typename idx_array_type>
inline void kk_write_2Dview_to_file(idx_array_type view, const char *filename) {
  typedef typename idx_array_type::HostMirror host_type;
  // typedef typename idx_array_type::size_type idx;
  host_type host_view = Kokkos::create_mirror_view(view);
  Kokkos::deep_copy(host_view, view);
  Kokkos::fence();
  std::ofstream myFile(filename, std::ios::out);
  for (size_t i = 0; i < view.extent(0); ++i) {
    for (size_t j = 0; j < view.extent(1); ++j) {
      myFile << host_view(i, j) << " ";
    }
    myFile << std::endl;
  }
  myFile.close();
}

template <typename idx_array_type>
inline void kk_read_2Dview_from_file(idx_array_type &view, const char *filename) {
  typedef typename idx_array_type::HostMirror host_type;
  // typedef typename idx_array_type::size_type idx;
  host_type host_view = Kokkos::create_mirror_view(view);
  std::ifstream myFile(filename, std::ios::in);

  for (size_t i = 0; i < view.extent(0); ++i) {
    for (size_t j = 0; j < view.extent(1); ++j) {
      myFile >> host_view(i, j);
    }
  }
  myFile.close();
  Kokkos::deep_copy(view, host_view);
  Kokkos::fence();
}

template <typename idx_array_type>
inline void kk_write_3Dview_to_file(idx_array_type view, const char *filename) {
  typedef typename idx_array_type::HostMirror host_type;
  // typedef typename idx_array_type::size_type idx;
  host_type host_view = Kokkos::create_mirror_view(view);
  Kokkos::deep_copy(host_view, view);
  Kokkos::fence();
  std::ofstream myFile(filename, std::ios::out);
  for (size_t i = 0; i < view.extent(0); ++i) {
    for (size_t j = 0; j < view.extent(1); ++j) {
      for (size_t k = 0; k < view.extent(2); ++k) {
        myFile << host_view(i, j, k) << " ";
      }
      myFile << std::endl;
    }
    myFile << std::endl;
  }
  myFile.close();
}

template <typename idx_array_type>
inline void kk_read_3Dview_from_file(idx_array_type &view, const char *filename) {
  typedef typename idx_array_type::HostMirror host_type;
  // typedef typename idx_array_type::size_type idx;
  host_type host_view = Kokkos::create_mirror_view(view);
  std::ifstream myFile(filename, std::ios::in);

  for (size_t i = 0; i < view.extent(0); ++i) {
    for (size_t j = 0; j < view.extent(1); ++j) {
      for (size_t k = 0; k < view.extent(2); ++k) {
        myFile >> host_view(i, j, k);
      }
    }
  }
  myFile.close();
  Kokkos::deep_copy(view, host_view);
  Kokkos::fence();
}

template <typename idx, typename wt>
[[deprecated]] void write_edgelist_bin(size_t ne, const idx *edge_begins, const idx *edge_ends, const wt *ew,
                                       const char *filename) {
  std::ofstream myFile(filename, std::ios::out | std::ios::binary);
  myFile.write((char *)&ne, sizeof(idx));
  myFile.write((char *)edge_begins, sizeof(idx) * (ne));
  myFile.write((char *)edge_ends, sizeof(idx) * (ne));
  myFile.write((char *)ew, sizeof(wt) * (ne));
  myFile.close();
}

template <typename idx, typename wt>
void read_edgelist_bin(idx *ne, idx **edge_begins, idx **edge_ends, wt **ew, const char *filename) {
  std::ifstream myFile(filename, std::ios::in | std::ios::binary);

  myFile.read((char *)ne, sizeof(idx));
  md_malloc<idx>(edge_begins, *ne);
  md_malloc<idx>(edge_ends, *ne);
  md_malloc<wt>(ew, *ne);
  myFile.read((char *)*edge_begins, sizeof(idx) * (*ne));
  myFile.read((char *)*edge_ends, sizeof(idx) * (*ne));
  myFile.read((char *)*ew, sizeof(wt) * (*ne));
  myFile.close();
}

inline bool endswith(std::string const &fullString, std::string const &ending) {
  if (fullString.length() >= ending.length()) {
    return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
  } else {
    return false;
  }
}

}  // namespace Impl
}  // namespace KokkosKernels

#endif
