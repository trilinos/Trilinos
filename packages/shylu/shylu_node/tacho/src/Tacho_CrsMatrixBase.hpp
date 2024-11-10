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
#ifndef __TACHO_CRS_MATRIX_BASE_HPP__
#define __TACHO_CRS_MATRIX_BASE_HPP__

/// \file Tacho_CrsMatrixBase.hpp
/// \brief CRS matrix base object interfaces to user provided input matrices.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho.hpp"

namespace Tacho {

/// \class CrsMatrixBase
/// \breif CRS matrix base object using Kokkos view and subview
template <typename ValueType, typename DeviceType> class CrsMatrixBase {
public:
  typedef ValueType value_type;
  typedef DeviceType device_type;

  typedef typename device_type::execution_space exec_space;
  typedef typename device_type::memory_space exec_memory_space;

  typedef Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace> host_device_type;

  typedef typename host_device_type::execution_space host_space;
  typedef typename host_device_type::memory_space host_memory_space;

  typedef Kokkos::View<value_type *, device_type> value_type_array;
  typedef Kokkos::View<ordinal_type *, device_type> ordinal_type_array;
  typedef Kokkos::View<size_type *, device_type> size_type_array;

  template <typename, typename> friend class CrsMatrixBase;

private:
  ordinal_type _m; //!< # of rows
  ordinal_type _n; //!< # of cols
  size_type _nnz;  //!< # of nonzeros

  size_type_array _ap;    //!< pointers to column index and values
  ordinal_type_array _aj; //!< column index compressed format
  value_type_array _ax;   //!< values

protected:
  void createInternal(const ordinal_type m, const ordinal_type n, const size_type nnz) {
    _m = m;
    _n = n;
    _nnz = nnz;

    if (static_cast<ordinal_type>(_ap.extent(0)) < (m + 1))
      _ap = size_type_array("CrsMatrixBase::RowPtrArray", m + 1);
    else
      Kokkos::deep_copy(_ap, size_type());

    if (static_cast<size_type>(_aj.extent(0)) < nnz)
      _aj = ordinal_type_array("CrsMatrixBase::ColsArray", nnz);
    else
      Kokkos::deep_copy(_aj, ordinal_type());

    if (static_cast<size_type>(_ax.extent(0)) < nnz)
      _ax = value_type_array("CrsMatrixBase::ValuesArray", nnz);
    else
      Kokkos::deep_copy(_ax, value_type());
  }

public:
  /// Interface functions
  /// ------------------------------------------------------------------
  void setExternalMatrix(const ordinal_type m, const ordinal_type n, const size_type nnz, const size_type_array &ap,
                         const ordinal_type_array &aj, const value_type_array &ax) {
    _m = m;
    _n = n;
    _nnz = nnz;
    _ap = ap;
    _aj = aj;
    _ax = ax;
  }

  void setExternalMatrix(const ordinal_type m, const ordinal_type n, const ordinal_type nnz, const size_type *ap,
                         const ordinal_type *aj, const value_type *ax) {
    _m = m;
    _n = n;
    _nnz = nnz;
    _ap = size_type_array(ap, _m + 1);
    _aj = ordinal_type_array(aj, _nnz);
    _ax = value_type_array(ax, _nnz);
  }

  void setNumNonZeros() {
    if (_m) {
      auto last = Kokkos::subview(_ap, _m);
      auto h_last = Kokkos::create_mirror_view(host_memory_space(), last);
      Kokkos::deep_copy(h_last, last);
      _nnz = h_last();
    }
  }

  KOKKOS_INLINE_FUNCTION
  size_type_array &RowPtr() { return _ap; }

  KOKKOS_INLINE_FUNCTION
  ordinal_type_array &Cols() { return _aj; }

  KOKKOS_INLINE_FUNCTION
  value_type_array &Values() { return _ax; }

  KOKKOS_INLINE_FUNCTION
  ordinal_type NumRows() const { return _m; }

  KOKKOS_INLINE_FUNCTION
  ordinal_type NumCols() const { return _n; }

  KOKKOS_INLINE_FUNCTION
  size_type NumNonZeros() const { return _nnz; }

  KOKKOS_INLINE_FUNCTION
  size_type &RowPtrBegin(const ordinal_type i) { return _ap[i]; }

  KOKKOS_INLINE_FUNCTION
  size_type RowPtrBegin(const ordinal_type i) const { return _ap[i]; }

  KOKKOS_INLINE_FUNCTION
  size_type &RowPtrEnd(const ordinal_type i) { return _ap[i + 1]; }

  KOKKOS_INLINE_FUNCTION
  size_type RowPtrEnd(const ordinal_type i) const { return _ap[i + 1]; }

  KOKKOS_INLINE_FUNCTION
  ordinal_type &Col(const ordinal_type k) { return _aj[k]; }

  KOKKOS_INLINE_FUNCTION
  ordinal_type Col(const ordinal_type k) const { return _aj[k]; }

  KOKKOS_INLINE_FUNCTION
  value_type &Value(const ordinal_type k) { return _ax[k]; }

  KOKKOS_INLINE_FUNCTION
  value_type Value(const ordinal_type k) const { return _ax[k]; }

  /// Constructors
  /// ------------------------------------------------------------------

  /// \brief Default constructor.
  CrsMatrixBase() : _m(0), _n(0), _nnz(0), _ap(), _aj(), _ax() {}

  /// \brief Constructor with label
  CrsMatrixBase(const ordinal_type m, const ordinal_type n, const size_type nnz)
      : _m(0), _n(0), _nnz(0), _ap(), _aj(), _ax() {
    createInternal(m, n, nnz);
  }

  /// \brief Copy constructor (shallow copy), for deep-copy use a method copy
  KOKKOS_INLINE_FUNCTION
  CrsMatrixBase(const CrsMatrixBase &b) : _m(b._m), _n(b._n), _nnz(b._nnz), _ap(b._ap), _aj(b._aj), _ax(b._ax) {}

  KOKKOS_INLINE_FUNCTION
  ~CrsMatrixBase() {} //= default;

  /// Create
  /// ------------------------------------------------------------------

  void clear() {
    _m = 0;
    _n = 0;
    _nnz = 0;

    _ap = size_type_array();
    _aj = ordinal_type_array();
    _ax = value_type_array();
  }

  void create(const ordinal_type m, const ordinal_type n, const size_type nnz) { createInternal(m, n, nnz); }

  template <typename SpT> void createConfTo(const CrsMatrixBase<value_type, SpT> &b) {
    createInternal(b._m, b._n, b._nnz);
  }

  /// Create
  /// ------------------------------------------------------------------

  /// \brief deep copy of matrix b
  template <typename SpT> inline void createMirror(const CrsMatrixBase<value_type, SpT> &b) {
    _m = b._m;
    _n = b._n;
    _nnz = b._nnz;

    _ap = Kokkos::create_mirror_view(exec_memory_space(), b._ap);
    _aj = Kokkos::create_mirror_view(exec_memory_space(), b._aj);
    _ax = Kokkos::create_mirror_view(exec_memory_space(), b._ax);
  }

  /// \brief deep copy of matrix b
  template <typename SpT> inline void copy(const CrsMatrixBase<value_type, SpT> &b) {
    Kokkos::deep_copy(_ap, b._ap);
    Kokkos::deep_copy(_aj, b._aj);
    Kokkos::deep_copy(_ax, b._ax);
  }

  /// \brief print out to stream
  inline std::ostream &showMe(std::ostream &os, const bool detail = false) const {
    std::streamsize prec = os.precision();
    os.precision(16);
    os << std::scientific;

    os << " -- CrsMatrixBase -- " << std::endl
       << "    # of Rows          = " << _m << std::endl
       << "    # of Cols          = " << _n << std::endl
       << "    # of NonZeros      = " << _nnz << std::endl
       << std::endl
       << "    RowPtrArray length = " << _ap.extent(0) << std::endl
       << "    ColArray length    = " << _aj.extent(0) << std::endl
       << "    ValueArray length  = " << _ax.extent(0) << std::endl
       << std::endl
       << "    Memory             = "
       << double(_ap.span() * sizeof(size_type) + _aj.span() * sizeof(ordinal_type) + _ax.span() * sizeof(value_type)) /
              1e6
       << " MB" << std::endl
       << std::endl;

    if (detail) {
      const int w = 10;
      if ((ordinal_type(_ap.size()) > _m) && (size_type(_aj.size()) >= _nnz) && (size_type(_ax.size()) >= _nnz)) {
        os << std::setw(w) << "Row"
           << "  " << std::setw(w) << "Col"
           << "  " << std::setw(w) << "Val" << std::endl;
        auto h_ap = Kokkos::create_mirror_view_and_copy(host_memory_space(), _ap);
        auto h_aj = Kokkos::create_mirror_view_and_copy(host_memory_space(), _aj);
        auto h_ax = Kokkos::create_mirror_view_and_copy(host_memory_space(), _ax);
        for (ordinal_type i = 0; i < _m; ++i) {
          const size_type jbegin = h_ap[i], jend = h_ap[i + 1];
          for (size_type j = jbegin; j < jend; ++j) {
            value_type val = h_ax[j];
            os << std::setw(w) << i << "  " << std::setw(w) << h_aj[j] << "  " << std::setw(w) << std::showpos << val
               << std::noshowpos << std::endl;
          }
        }
      }
    }

    os.unsetf(std::ios::scientific);
    os.precision(prec);

    return os;
  }
};

// A = P B P^{-1}
template <typename CrsMatrixType, typename OrdinalTypeArray>
inline static void applyPermutationToCrsMatrix(/* */ CrsMatrixType &A, const CrsMatrixType &B,
                                               const OrdinalTypeArray &p, const OrdinalTypeArray &ip) {
  const ordinal_type m = A.NumRows(); //, n = A.NumCols();
  typedef typename CrsMatrixType::exec_space exec_space;
  typedef typename CrsMatrixType::exec_memory_space exec_memory_space;

  /// temporary matrix with the same structure
  auto ap = A.RowPtr();
  auto aj = A.Cols();
  auto ax = A.Values();

  auto perm = Kokkos::create_mirror_view(exec_memory_space(), p);
  Kokkos::deep_copy(perm, p);
  auto peri = Kokkos::create_mirror_view(exec_memory_space(), ip);
  Kokkos::deep_copy(peri, ip);
  { /// permute row indices (exclusive scan)
    Kokkos::RangePolicy<exec_space, Kokkos::Schedule<Kokkos::Static>> policy(0, m + 1);
    Kokkos::parallel_scan(
        policy, KOKKOS_LAMBDA(const ordinal_type &i, size_type &update, const bool &final) {
          if (final)
            ap(i) = update;

          if (i < m) {
            const ordinal_type row = perm(i);
            update += (B.RowPtrEnd(row) - B.RowPtrBegin(row));
          }
        });
    Kokkos::fence();
  }
  { /// permute col indices (do not sort)
    typedef Kokkos::TeamPolicy<exec_space> team_policy_type;
    team_policy_type policy(m, 1, 1); // Kokkos::AUTO()); ///, Kokkos::AUTO());
    Kokkos::parallel_for(
        policy, KOKKOS_LAMBDA(const typename team_policy_type::member_type &member) {
          const ordinal_type i = member.league_rank(), /// row in A
              kbeg = ap(i), kend = ap(i + 1), nk = kend - kbeg,
                             row = perm(i), /// row in B
              colbeg = B.RowPtrBegin(row);
          // colend = B.RowPtrEnd(row);

          Kokkos::parallel_for(
              Kokkos::TeamVectorRange(member, nk),
              [&, kbeg, colbeg](const ordinal_type &k) { /// compiler bug with c++14 lambda capturing and workaround
                const ordinal_type tk = kbeg + k, sk = colbeg + k;
                aj(tk) = peri(B.Col(sk));
                ax(tk) = B.Value(sk);
              });
        });
    Kokkos::fence();
  }
}

// A = P B P^{-1}
template <typename CrsMatrixType, typename OrdinalTypeArray>
inline static void applyPermutationToCrsMatrixLower(/* */ CrsMatrixType &A, const CrsMatrixType &B,
                                                    const OrdinalTypeArray &p, const OrdinalTypeArray &ip) {
  const ordinal_type m = A.NumRows(); //, n = A.NumCols();
  typedef typename CrsMatrixType::exec_space exec_space;
  typedef typename CrsMatrixType::exec_memory_space exec_memory_space;

  /// temporary matrix with the same structure
  auto ap = A.RowPtr();
  auto aj = A.Cols();
  auto ax = A.Values();

  auto perm = Kokkos::create_mirror_view(exec_memory_space(), p);
  Kokkos::deep_copy(perm, p);
  auto peri = Kokkos::create_mirror_view(exec_memory_space(), ip);
  Kokkos::deep_copy(peri, ip);

  { /// permute row indices (exclusive scan)
    Kokkos::RangePolicy<exec_space, Kokkos::Schedule<Kokkos::Static>> policy(0, m + 1);
    Kokkos::parallel_scan(
        policy, KOKKOS_LAMBDA(const ordinal_type &i, size_type &update, const bool &final) {
          if (final)
            ap(i) = update;

          if (i < m) {
            ordinal_type count(0);
            const ordinal_type row = perm(i); /// row in B
            for (ordinal_type k = B.RowPtrBegin(row), kend = B.RowPtrEnd(row); k < kend; ++k) {
              const ordinal_type j = peri(B.Col(k)); /// col in A
              count += (i >= j);                     /// lower triangular
            }
            update += count;
          }
        });
    Kokkos::fence();
  }
  { /// permute col indices (do not sort)
    Kokkos::RangePolicy<exec_space, Kokkos::Schedule<Kokkos::Static>> policy(0, m);
    Kokkos::parallel_for(
        policy, KOKKOS_LAMBDA(const ordinal_type &i) { /// row in A
          const ordinal_type kbeg = ap(i),
                             row = perm(i), /// row in B
              colbeg = B.RowPtrBegin(row), colend = B.RowPtrEnd(row), nk = colend - colbeg;

          for (ordinal_type k = 0, t = 0; k < nk; ++k) {
            const ordinal_type tk = kbeg + t, sk = colbeg + k;
            const ordinal_type j = peri(B.Col(sk)); /// col in A
            if (i >= j) {
              aj(tk) = j;
              ax(tk) = B.Value(sk);
              ++t;
            }
          }
        });
    Kokkos::fence();
  }

  /// update A with new nnz (for now let's not resize)
  A.setNumNonZeros();
}

template <typename ValueType, typename DeviceType>
inline double computeRelativeResidual(const CrsMatrixBase<ValueType, DeviceType> &A,
                                      const Kokkos::View<ValueType **, Kokkos::LayoutLeft, DeviceType> &x,
                                      const Kokkos::View<ValueType **, Kokkos::LayoutLeft, DeviceType> &b) {
  const bool test = (size_t(A.NumRows()) != size_t(A.NumCols()) || size_t(A.NumRows()) != size_t(b.extent(0)) ||
                     size_t(x.extent(0)) != size_t(b.extent(0)) || size_t(x.extent(1)) != size_t(b.extent(1)));
  if (test)
    throw std::logic_error("A,x and b dimensions are not compatible");

  typedef ValueType value_type;
  typedef typename UseThisDevice<Kokkos::DefaultHostExecutionSpace>::type host_device_type;
  typedef typename host_device_type::memory_space host_memory_space;

  CrsMatrixBase<value_type, host_device_type> h_A;
  h_A.createMirror(A);
  h_A.copy(A);

  auto h_x = Kokkos::create_mirror_view(host_memory_space(), x);
  Kokkos::deep_copy(h_x, x);
  auto h_b = Kokkos::create_mirror_view(host_memory_space(), b);
  Kokkos::deep_copy(h_b, b);

  typedef ArithTraits<value_type> arith_traits;
  const ordinal_type m = h_A.NumRows(), k = h_b.extent(1);
  double diff = 0, norm = 0;
  for (ordinal_type i = 0; i < m; ++i) {
    for (ordinal_type p = 0; p < k; ++p) {
      value_type s = 0;
      const ordinal_type jbeg = h_A.RowPtrBegin(i), jend = h_A.RowPtrEnd(i);
      for (ordinal_type j = jbeg; j < jend; ++j) {
        const ordinal_type col = h_A.Col(j);
        s += h_A.Value(j) * h_x(col, p);
      }
      norm += arith_traits::real(h_b(i, p) * arith_traits::conj(h_b(i, p)));
      diff += arith_traits::real((h_b(i, p) - s) * arith_traits::conj(h_b(i, p) - s));
    }
  }
  return sqrt(diff / norm);
}

template <typename ValueType, typename DeviceType>
inline void computeSpMV(const CrsMatrixBase<ValueType, DeviceType> &A,
                        const Kokkos::View<ValueType **, Kokkos::LayoutLeft, DeviceType> &x,
                              Kokkos::View<ValueType **, Kokkos::LayoutLeft, DeviceType> &b) {
  const bool test = (size_t(A.NumRows()) != size_t(A.NumCols()) || size_t(A.NumRows()) != size_t(b.extent(0)) ||
                     size_t(x.extent(0)) != size_t(b.extent(0)) || size_t(x.extent(1)) != size_t(b.extent(1)));
  if (test)
    throw std::logic_error("A,x and b dimensions are not compatible");

  typedef ValueType value_type;
  typedef typename UseThisDevice<Kokkos::DefaultHostExecutionSpace>::type host_device_type;
  typedef typename host_device_type::memory_space host_memory_space;

  CrsMatrixBase<value_type, host_device_type> h_A;
  h_A.createMirror(A);
  h_A.copy(A);

  auto h_x = Kokkos::create_mirror_view(host_memory_space(), x);
  Kokkos::deep_copy(h_x, x);
  auto h_b = Kokkos::create_mirror_view(host_memory_space(), b);

  const ordinal_type m = h_A.NumRows(), k = h_b.extent(1);
  for (ordinal_type i = 0; i < m; ++i) {
    for (ordinal_type p = 0; p < k; ++p) {
      value_type s = 0;
      const ordinal_type jbeg = h_A.RowPtrBegin(i), jend = h_A.RowPtrEnd(i);
      for (ordinal_type j = jbeg; j < jend; ++j) {
        const ordinal_type col = h_A.Col(j);
        s += h_A.Value(j) * h_x(col, p);
      }
      h_b(i,p) = s;
    }
  }
  Kokkos::deep_copy(b, h_b);
  return;
}

} // namespace Tacho
#endif
