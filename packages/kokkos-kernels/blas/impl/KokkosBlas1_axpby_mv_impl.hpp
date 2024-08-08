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
#ifndef KOKKOSBLAS1_AXPBY_MV_IMPL_HPP_
#define KOKKOSBLAS1_AXPBY_MV_IMPL_HPP_

#include <KokkosBlas1_axpby_impl.hpp>

namespace KokkosBlas {
namespace Impl {

//
// axpby
//

// Functor for multivectors X and Y and 1-D views a and b, that
// computes any of the following:
//
// 1. Y(i,j) = alpha*X(i,j) + beta*Y(i,j) for alpha,beta in -1,0,1
// 2. Y(i,j) = a(j)*X(i,j) + beta*Y(i,j) for beta in -1,0,1
// 3. Y(i,j) = alpha*X(i,j) + beta*Y(i,j) for alpha in -1,0,1
// 4. Y(i,j) = a(j)*X(i,j) + b(j)*Y(i,j)
//
// The template parameters scalar_x and scalar_y correspond to alpha
// resp. beta in the operation y = alpha*x + beta*y.  The values -1,
// 0, and -1 correspond to literal values of those coefficients.
// The value 2 tells the functor to use the corresponding vector of
// coefficients.  Any literal coefficient of zero has BLAS semantics
// of ignoring the corresponding (multi)vector entry.  This does not
// apply to coefficients in the a and b vectors, if they are used.
template <class AV, class XMV, class BV, class YMV, int scalar_x, int scalar_y,
          class SizeType = typename YMV::size_type>
struct Axpby_MV_Functor {
  typedef SizeType size_type;
  typedef Kokkos::ArithTraits<typename YMV::non_const_value_type> ATS;

  const size_type numCols;
  XMV m_x;
  YMV m_y;
  AV m_a;
  BV m_b;

  Axpby_MV_Functor(const XMV& X, const YMV& Y, const AV& av, const BV& bv)
      : numCols(X.extent(1)), m_x(X), m_y(Y), m_a(av), m_b(bv) {
    static_assert(Kokkos::is_view<AV>::value,
                  "KokkosBlas::Impl::Axpby_MV_Functor(ABgeneric)"
                  ": 'a' is not a Kokkos::View.");
    static_assert(Kokkos::is_view<XMV>::value,
                  "KokkosBlas::Impl::Axpby_MV_Functor(ABgeneric)"
                  ": X is not a Kokkos::View.");
    static_assert(Kokkos::is_view<BV>::value,
                  "KokkosBlas::Impl::Axpby_MV_Functor(ABgeneric)"
                  ": 'b' is not a Kokkos::View.");
    static_assert(Kokkos::is_view<YMV>::value,
                  "KokkosBlas::Impl::Axpby_MV_Functor(ABgeneric)"
                  ": Y is not a Kokkos::View.");
    static_assert(std::is_same<typename YMV::value_type, typename YMV::non_const_value_type>::value,
                  "KokkosBlas::Impl::Axpby_MV_Functor(ABgeneric)"
                  ": Y must be nonconst, since it is an output argument"
                  " and we have to be able to write to its entries.");
    static_assert((int)YMV::rank == (int)XMV::rank,
                  "KokkosBlas::Impl::Axpby_MV_Functor(ABgeneric)"
                  ": X and Y must have the same rank.");
    static_assert(YMV::rank == 2,
                  "KokkosBlas::Impl::Axpby_MV_Functor(ABgeneric)"
                  ": XMV and YMV must have rank 2.");
    static_assert(AV::rank == 1,
                  "KokkosBlas::Impl::Axpby_MV_Functor(ABgeneric)"
                  ": AV must have rank 1.");
    static_assert(BV::rank == 1,
                  "KokkosBlas::Impl::Axpby_MV_Functor(ABgeneric)"
                  ": BV must have rank 1.");
    static_assert((-1 <= scalar_x) && (scalar_x <= 2) && (-1 <= scalar_y) && (scalar_y <= 2),
                  "KokkosBlas::Impl::Axpby_MV_Functor(ABgeneric)"
                  ": scalar_x and/or scalar_y are out of range.");
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_type& i) const {
    // scalar_x and scalar_y are compile-time constants (since they
    // are template parameters), so the compiler should evaluate these
    // branches at compile time.

    // **************************************************************
    // Possibilities with 'scalar_x == 0'
    // **************************************************************
    if constexpr (scalar_x == 0) {
      if constexpr (scalar_y == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
        for (size_type k = 0; k < numCols; ++k) {
          m_y(i, k) = ATS::zero();
        }
      } else if constexpr (scalar_y == -1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
        for (size_type k = 0; k < numCols; ++k) {
          m_y(i, k) = -m_y(i, k);
        }
      } else if constexpr (scalar_y == 1) {
        // Nothing to do: Y(i,j) := Y(i,j)
      } else if constexpr (scalar_y == 2) {
        if (m_b.extent(0) == 1) {
          if (m_b(0) == Kokkos::ArithTraits<typename BV::non_const_value_type>::zero()) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
            for (size_type k = 0; k < numCols; ++k) {
              m_y(i, k) = Kokkos::ArithTraits<typename YMV::non_const_value_type>::zero();
            }
          } else {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
            for (size_type k = 0; k < numCols; ++k) {
              m_y(i, k) = m_b(0) * m_y(i, k);
            }
          }
        } else {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
          for (size_type k = 0; k < numCols; ++k) {
            m_y(i, k) = m_b(k) * m_y(i, k);
          }
        }
      }
    }
    // **************************************************************
    // Possibilities with 'scalar_x == -1'
    // **************************************************************
    else if constexpr (scalar_x == -1) {
      if constexpr (scalar_y == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
        for (size_type k = 0; k < numCols; ++k) {
          m_y(i, k) = -m_x(i, k);
        }
      } else if constexpr (scalar_y == -1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
        for (size_type k = 0; k < numCols; ++k) {
          m_y(i, k) = -m_x(i, k) - m_y(i, k);
        }
      } else if constexpr (scalar_y == 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
        for (size_type k = 0; k < numCols; ++k) {
          m_y(i, k) = -m_x(i, k) + m_y(i, k);
        }
      } else if constexpr (scalar_y == 2) {
        if (m_b.extent(0) == 1) {
          if (m_b(0) == Kokkos::ArithTraits<typename BV::non_const_value_type>::zero()) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
            for (size_type k = 0; k < numCols; ++k) {
              m_y(i, k) = -m_x(i, k);
            }
          } else {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
            for (size_type k = 0; k < numCols; ++k) {
              m_y(i, k) = -m_x(i, k) + m_b(0) * m_y(i, k);
            }
          }
        } else {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
          for (size_type k = 0; k < numCols; ++k) {
            m_y(i, k) = -m_x(i, k) + m_b(k) * m_y(i, k);
          }
        }
      }
    }
    // **************************************************************
    // Possibilities with 'scalar_x == 1'
    // **************************************************************
    else if constexpr (scalar_x == 1) {
      if constexpr (scalar_y == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
        for (size_type k = 0; k < numCols; ++k) {
          m_y(i, k) = m_x(i, k);
        }
      } else if constexpr (scalar_y == -1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
        for (size_type k = 0; k < numCols; ++k) {
          m_y(i, k) = m_x(i, k) - m_y(i, k);
        }
      } else if constexpr (scalar_y == 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
        for (size_type k = 0; k < numCols; ++k) {
          m_y(i, k) = m_x(i, k) + m_y(i, k);
        }
      } else if constexpr (scalar_y == 2) {
        if (m_b.extent(0) == 1) {
          if (m_b(0) == Kokkos::ArithTraits<typename BV::non_const_value_type>::zero()) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
            for (size_type k = 0; k < numCols; ++k) {
              m_y(i, k) = m_x(i, k);
            }
          } else {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
            for (size_type k = 0; k < numCols; ++k) {
              m_y(i, k) = m_x(i, k) + m_b(0) * m_y(i, k);
            }
          }
        } else {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
          for (size_type k = 0; k < numCols; ++k) {
            m_y(i, k) = m_x(i, k) + m_b(k) * m_y(i, k);
          }
        }
      }
    }
    // **************************************************************
    // Possibilities with 'scalar_x == 2'
    // **************************************************************
    else if constexpr (scalar_x == 2) {
      if constexpr (scalar_y == 0) {
        if (m_a.extent(0) == 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
          for (size_type k = 0; k < numCols; ++k) {
            m_y(i, k) = m_a(0) * m_x(i, k);
          }
        } else {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
          for (size_type k = 0; k < numCols; ++k) {
            m_y(i, k) = m_a(k) * m_x(i, k);
          }
        }
      } else if constexpr (scalar_y == -1) {
        if (m_a.extent(0) == 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
          for (size_type k = 0; k < numCols; ++k) {
            m_y(i, k) = m_a(0) * m_x(i, k) - m_y(i, k);
          }
        } else {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
          for (size_type k = 0; k < numCols; ++k) {
            m_y(i, k) = m_a(k) * m_x(i, k) - m_y(i, k);
          }
        }
      } else if constexpr (scalar_y == 1) {
        if (m_a.extent(0) == 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
          for (size_type k = 0; k < numCols; ++k) {
            m_y(i, k) = m_a(0) * m_x(i, k) + m_y(i, k);
          }
        } else {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
          for (size_type k = 0; k < numCols; ++k) {
            m_y(i, k) = m_a(k) * m_x(i, k) + m_y(i, k);
          }
        }
      } else if constexpr (scalar_y == 2) {
        if (m_a.extent(0) == 1) {
          if (m_b.extent(0) == 1) {
            if (m_b(0) == Kokkos::ArithTraits<typename BV::non_const_value_type>::zero()) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
              for (size_type k = 0; k < numCols; ++k) {
                m_y(i, k) = m_a(0) * m_x(i, k);
              }
            } else {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
              for (size_type k = 0; k < numCols; ++k) {
                m_y(i, k) = m_a(0) * m_x(i, k) + m_b(0) * m_y(i, k);
              }
            }
          } else {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
            for (size_type k = 0; k < numCols; ++k) {
              m_y(i, k) = m_a(0) * m_x(i, k) + m_b(k) * m_y(i, k);
            }
          }
        } else {
          if (m_b.extent(0) == 1) {
            if (m_b(0) == Kokkos::ArithTraits<typename BV::non_const_value_type>::zero()) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
              for (size_type k = 0; k < numCols; ++k) {
                m_y(i, k) = m_a(k) * m_x(i, k);
              }
            } else {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
              for (size_type k = 0; k < numCols; ++k) {
                m_y(i, k) = m_a(k) * m_x(i, k) + m_b(0) * m_y(i, k);
              }
            }
          } else {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
            for (size_type k = 0; k < numCols; ++k) {
              m_y(i, k) = m_a(k) * m_x(i, k) + m_b(k) * m_y(i, k);
            }
          }
        }
      }  // if constexpr (scalar_y == ...) else if
    }    // if constexpr (scalar_x == ...) else if
  }      // void operator()
};

// Variant of Axpby_MV_Functor, where a and b are scalars.
// This functor computes any of the following:
//
// 1. Y(i,j) = alpha*X(i,j) + beta*Y(i,j) for alpha,beta in -1,0,1
// 2. Y(i,j) = a*X(i,j) + beta*Y(i,j) for beta in -1,0,1
// 3. Y(i,j) = alpha*X(i,j) + b*Y(i,j) for alpha in -1,0,1
// 4. Y(i,j) = a*X(i,j) + b*Y(i,j)
//
// The template parameters scalar_x and scalar_y correspond to alpha
// resp. beta in the operation y = alpha*x + beta*y.  The values -1,
// 0, and -1 correspond to literal values of those coefficients.
// The value 2 tells the functor to use the corresponding vector of
// coefficients.  Any literal coefficient of zero has BLAS semantics
// of ignoring the corresponding (multi)vector entry.  This does not
// apply to coefficients in the a and b vectors, if they are used.
//
// This version works by partial specialization on AV and BV.
// In this partial specialization, both AV and BV are scalars.
template <class XMV, class YMV, int scalar_x, int scalar_y, class SizeType>
struct Axpby_MV_Functor<typename XMV::non_const_value_type, XMV, typename YMV::non_const_value_type, YMV, scalar_x,
                        scalar_y, SizeType> {
  typedef SizeType size_type;
  typedef Kokkos::ArithTraits<typename YMV::non_const_value_type> ATS;

  const size_type numCols;
  XMV m_x;
  YMV m_y;
  const typename XMV::non_const_value_type m_a;
  const typename YMV::non_const_value_type m_b;

  Axpby_MV_Functor(const XMV& X, const YMV& Y, const typename XMV::non_const_value_type& a,
                   const typename YMV::non_const_value_type& b)
      : numCols(X.extent(1)), m_x(X), m_y(Y), m_a(a), m_b(b) {
    static_assert(Kokkos::is_view<XMV>::value,
                  "KokkosBlas::Impl::Axpby_MV_Functor(ABscalars)"
                  ": X is not a Kokkos::View.");
    static_assert(Kokkos::is_view<YMV>::value,
                  "KokkosBlas::Impl::Axpby_MV_Functor(ABscalars)"
                  ": Y is not a Kokkos::View.");
    static_assert(std::is_same<typename YMV::value_type, typename YMV::non_const_value_type>::value,
                  "KokkosBlas::Impl::Axpby_MV_Functor(ABscalars)"
                  ": Y must be nonconst, since it is an output argument"
                  " and we have to be able to write to its entries.");
    static_assert((int)YMV::rank == (int)XMV::rank,
                  "KokkosBlas::Impl::Axpby_MV_Functor(ABscalars)"
                  ": X and Y must have the same rank.");
    static_assert(YMV::rank == 2,
                  "KokkosBlas::Impl::Axpby_MV_Functor(ABscalars)"
                  ": XMV and YMV must have rank 2.");
    static_assert((-1 <= scalar_x) && (scalar_x <= 2) && (-1 <= scalar_y) && (scalar_y <= 2),
                  "KokkosBlas::Impl::Axpby_MV_Functor(ABscalars)"
                  ": scalar_x and/or scalar_y are out of range.");
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_type& i) const {
    // scalar_x and scalar_y are compile-time constants (since they
    // are template parameters), so the compiler should evaluate these
    // branches at compile time.

    // **************************************************************
    // Possibilities with 'scalar_x == 0'
    // **************************************************************
    if constexpr (scalar_x == 0) {
      if constexpr (scalar_y == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
        for (size_type k = 0; k < numCols; ++k) {
          m_y(i, k) = ATS::zero();
        }
      } else if constexpr (scalar_y == -1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
        for (size_type k = 0; k < numCols; ++k) {
          m_y(i, k) = -m_y(i, k);
        }
      } else if constexpr (scalar_y == 1) {
        // Nothing to do: Y(i,j) := Y(i,j)
      } else if constexpr (scalar_y == 2) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
        for (size_type k = 0; k < numCols; ++k) {
          m_y(i, k) = m_b * m_y(i, k);
        }
      }
    }
    // **************************************************************
    // Possibilities with 'scalar_x == -1'
    // **************************************************************
    else if constexpr (scalar_x == -1) {
      if constexpr (scalar_y == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
        for (size_type k = 0; k < numCols; ++k) {
          m_y(i, k) = -m_x(i, k);
        }
      } else if constexpr (scalar_y == -1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
        for (size_type k = 0; k < numCols; ++k) {
          m_y(i, k) = -m_x(i, k) - m_y(i, k);
        }
      } else if constexpr (scalar_y == 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
        for (size_type k = 0; k < numCols; ++k) {
          m_y(i, k) = -m_x(i, k) + m_y(i, k);
        }
      } else if constexpr (scalar_y == 2) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
        for (size_type k = 0; k < numCols; ++k) {
          m_y(i, k) = -m_x(i, k) + m_b * m_y(i, k);
        }
      }
    }
    // **************************************************************
    // Possibilities with 'scalar_x == 1'
    // **************************************************************
    else if constexpr (scalar_x == 1) {
      if constexpr (scalar_y == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
        for (size_type k = 0; k < numCols; ++k) {
          m_y(i, k) = m_x(i, k);
        }
      } else if constexpr (scalar_y == -1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
        for (size_type k = 0; k < numCols; ++k) {
          m_y(i, k) = m_x(i, k) - m_y(i, k);
        }
      } else if constexpr (scalar_y == 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
        for (size_type k = 0; k < numCols; ++k) {
          m_y(i, k) = m_x(i, k) + m_y(i, k);
        }
      } else if constexpr (scalar_y == 2) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
        for (size_type k = 0; k < numCols; ++k) {
          m_y(i, k) = m_x(i, k) + m_b * m_y(i, k);
        }
      }
    }
    // **************************************************************
    // Possibilities with 'scalar_x == 2'
    // **************************************************************
    else if constexpr (scalar_x == 2) {
      if constexpr (scalar_y == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
        for (size_type k = 0; k < numCols; ++k) {
          m_y(i, k) = m_a * m_x(i, k);
        }
      } else if constexpr (scalar_y == -1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
        for (size_type k = 0; k < numCols; ++k) {
          m_y(i, k) = m_a * m_x(i, k) - m_y(i, k);
        }
      } else if constexpr (scalar_y == 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
        for (size_type k = 0; k < numCols; ++k) {
          m_y(i, k) = m_a * m_x(i, k) + m_y(i, k);
        }
      } else if constexpr (scalar_y == 2) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
        for (size_type k = 0; k < numCols; ++k) {
          m_y(i, k) = m_a * m_x(i, k) + m_b * m_y(i, k);
        }
      }  // if constexpr (scalar_y == ...) else if
    }    // if constexpr (scalar_x == ...) else if
  }      // void operator()
};

// Column-unrolled variant of Axpby_MV_Functor.  The number of columns
// in X and Y, UNROLL, is a compile-time constant.
template <class AV, class XMV, class BV, class YMV, int scalar_x, int scalar_y, int UNROLL, class SizeType>
struct Axpby_MV_Unroll_Functor {
  typedef SizeType size_type;
  typedef Kokkos::ArithTraits<typename YMV::non_const_value_type> ATS;

  XMV m_x;
  YMV m_y;
  AV m_a;
  BV m_b;

  Axpby_MV_Unroll_Functor(const XMV& x, const YMV& y, const AV& av, const BV& bv, const SizeType startingColumn)
      : m_x(x), m_y(y), m_a(av), m_b(bv) {
    static_assert(Kokkos::is_view<AV>::value,
                  "KokkosBlas::Impl::Axpby_MV_Unroll_Functor(ABgeneric)"
                  ": 'a' is not a Kokkos::View.");
    static_assert(Kokkos::is_view<XMV>::value,
                  "KokkosBlas::Impl::Axpby_MV_Unroll_Functor(ABgeneric)"
                  ": X is not a Kokkos::View.");
    static_assert(Kokkos::is_view<BV>::value,
                  "KokkosBlas::Impl::Axpby_MV_Unroll_Functor(ABgeneric)"
                  ": 'b' is not a Kokkos::View.");
    static_assert(Kokkos::is_view<YMV>::value,
                  "KokkosBlas::Impl::Axpby_MV_Unroll_Functor(ABgeneric)"
                  ": Y is not a Kokkos::View.");
    static_assert(std::is_same<typename YMV::value_type, typename YMV::non_const_value_type>::value,
                  "KokkosBlas::Impl::Axpby_MV_Unroll_Functor(ABgeneric)"
                  ": Y must be nonconst, since it is an output argument"
                  " and we have to be able to write to its entries.");
    static_assert((int)YMV::rank == (int)XMV::rank,
                  "KokkosBlas::Impl::Axpby_MV_Unroll_Functor(ABgeneric)"
                  ": X and Y must have the same rank.");
    static_assert(YMV::rank == 2,
                  "KokkosBlas::Impl::Axpby_MV_Unroll_Functor(ABgeneric)"
                  ": XMV and YMV must have rank 2.");
    static_assert(AV::rank == 1,
                  "KokkosBlas::Impl::Axpby_MV_Unroll_Functor(ABgeneric)"
                  ": AV must have rank 1.");
    static_assert(BV::rank == 1,
                  "KokkosBlas::Impl::Axpby_MV_Unroll_Functor(ABgeneric)"
                  ": BV must have rank 1.");
    static_assert((-1 <= scalar_x) && (scalar_x <= 2) && (-1 <= scalar_y) && (scalar_y <= 2),
                  "KokkosBlas::Impl::Axpby_MV_Unroll_Functor(ABgeneric)"
                  ": scalar_x and/or scalar_y are out of range.");

    if (startingColumn != 0) {
      if (axpbyVarExtent(m_a) > 1) {
        m_a = Kokkos::subview(av, std::make_pair(startingColumn, SizeType(av.extent(0))));
      }
      if (axpbyVarExtent(m_b) > 1) {
        m_b = Kokkos::subview(bv, std::make_pair(startingColumn, SizeType(bv.extent(0))));
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_type& i) const {
    // scalar_x and scalar_y are compile-time constants (since they
    // are template parameters), so the compiler should evaluate these
    // branches at compile time.

    // **************************************************************
    // Possibilities with 'scalar_x == 0'
    // **************************************************************
    if constexpr (scalar_x == 0) {
      if constexpr (scalar_y == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(i, k) = ATS::zero();
        }
      } else if constexpr (scalar_y == -1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(i, k) = -m_y(i, k);
        }
      } else if constexpr (scalar_y == 1) {
        // Nothing to do: Y(i,j) := Y(i,j)
      } else if constexpr (scalar_y == 2) {
        if (m_b.extent(0) == 1) {
          if (m_b(0) == Kokkos::ArithTraits<typename BV::non_const_value_type>::zero()) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
            for (int k = 0; k < UNROLL; ++k) {
              m_y(i, k) = Kokkos::ArithTraits<typename YMV::non_const_value_type>::zero();
            }
          } else {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
            for (int k = 0; k < UNROLL; ++k) {
              m_y(i, k) = m_b(0) * m_y(i, k);
            }
          }
        } else {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (int k = 0; k < UNROLL; ++k) {
            m_y(i, k) = m_b(k) * m_y(i, k);
          }
        }
      }
    }
    // **************************************************************
    // Possibilities with 'scalar_x == -1'
    // **************************************************************
    else if constexpr (scalar_x == -1) {
      if constexpr (scalar_y == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(i, k) = -m_x(i, k);
        }
      } else if constexpr (scalar_y == -1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(i, k) = -m_x(i, k) - m_y(i, k);
        }
      } else if constexpr (scalar_y == 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(i, k) = -m_x(i, k) + m_y(i, k);
        }
      } else if constexpr (scalar_y == 2) {
        if (m_b.extent(0) == 1) {
          if (m_b(0) == Kokkos::ArithTraits<typename BV::non_const_value_type>::zero()) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
            for (int k = 0; k < UNROLL; ++k) {
              m_y(i, k) = -m_x(i, k);
            }
          } else {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
            for (int k = 0; k < UNROLL; ++k) {
              m_y(i, k) = -m_x(i, k) + m_b(0) * m_y(i, k);
            }
          }
        } else {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (int k = 0; k < UNROLL; ++k) {
            m_y(i, k) = -m_x(i, k) + m_b(k) * m_y(i, k);
          }
        }
      }
    }
    // **************************************************************
    // Possibilities with 'scalar_x == 1'
    // **************************************************************
    else if constexpr (scalar_x == 1) {
      if constexpr (scalar_y == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(i, k) = m_x(i, k);
        }
      } else if constexpr (scalar_y == -1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(i, k) = m_x(i, k) - m_y(i, k);
        }
      } else if constexpr (scalar_y == 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(i, k) = m_x(i, k) + m_y(i, k);
        }
      } else if constexpr (scalar_y == 2) {
        if (m_b.extent(0) == 1) {
          if (m_b(0) == Kokkos::ArithTraits<typename BV::non_const_value_type>::zero()) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
            for (int k = 0; k < UNROLL; ++k) {
              m_y(i, k) = m_x(i, k);
            }
          } else {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
            for (int k = 0; k < UNROLL; ++k) {
              m_y(i, k) = m_x(i, k) + m_b(0) * m_y(i, k);
            }
          }
        } else {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (int k = 0; k < UNROLL; ++k) {
            m_y(i, k) = m_x(i, k) + m_b(k) * m_y(i, k);
          }
        }
      }
    }
    // **************************************************************
    // Possibilities with 'scalar_x == 2'
    // **************************************************************
    else if constexpr (scalar_x == 2) {
      if constexpr (scalar_y == 0) {
        if (m_a.extent(0) == 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (int k = 0; k < UNROLL; ++k) {
            m_y(i, k) = m_a(0) * m_x(i, k);
          }
        } else {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (int k = 0; k < UNROLL; ++k) {
            m_y(i, k) = m_a(k) * m_x(i, k);
          }
        }
      } else if constexpr (scalar_y == -1) {
        if (m_a.extent(0) == 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (int k = 0; k < UNROLL; ++k) {
            m_y(i, k) = m_a(0) * m_x(i, k) - m_y(i, k);
          }
        } else {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (int k = 0; k < UNROLL; ++k) {
            m_y(i, k) = m_a(k) * m_x(i, k) - m_y(i, k);
          }
        }
      } else if constexpr (scalar_y == 1) {
        if (m_a.extent(0) == 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (int k = 0; k < UNROLL; ++k) {
            m_y(i, k) = m_a(0) * m_x(i, k) + m_y(i, k);
          }
        } else {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (int k = 0; k < UNROLL; ++k) {
            m_y(i, k) = m_a(k) * m_x(i, k) + m_y(i, k);
          }
        }
      } else if constexpr (scalar_y == 2) {
        if (m_a.extent(0) == 1) {
          if (m_b.extent(0) == 1) {
            if (m_b(0) == Kokkos::ArithTraits<typename BV::non_const_value_type>::zero()) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
              for (int k = 0; k < UNROLL; ++k) {
                m_y(i, k) = m_a(0) * m_x(i, k);
              }
            } else {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
              for (int k = 0; k < UNROLL; ++k) {
                m_y(i, k) = m_a(0) * m_x(i, k) + m_b(0) * m_y(i, k);
              }
            }
          } else {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
            for (int k = 0; k < UNROLL; ++k) {
              m_y(i, k) = m_a(0) * m_x(i, k) + m_b(k) * m_y(i, k);
            }
          }
        } else {
          if (m_b.extent(0) == 1) {
            if (m_b(0) == Kokkos::ArithTraits<typename BV::non_const_value_type>::zero()) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
              for (int k = 0; k < UNROLL; ++k) {
                m_y(i, k) = m_a(k) * m_x(i, k);
              }
            } else {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
              for (int k = 0; k < UNROLL; ++k) {
                m_y(i, k) = m_a(k) * m_x(i, k) + m_b(0) * m_y(i, k);
              }
            }
          } else {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
            for (int k = 0; k < UNROLL; ++k) {
              m_y(i, k) = m_a(k) * m_x(i, k) + m_b(k) * m_y(i, k);
            }
          }
        }
      }
    }
  }
};

// Variant of Axpby_MV_Unroll_Functor for single coefficients (rather
// than vectors of coefficients) a and b.  The number of columns in X
// and Y, UNROLL, is a compile-time constant.
template <class XMV, class YMV, int scalar_x, int scalar_y, int UNROLL, class SizeType>
struct Axpby_MV_Unroll_Functor<typename XMV::non_const_value_type, XMV, typename YMV::non_const_value_type, YMV,
                               scalar_x, scalar_y, UNROLL, SizeType> {
  typedef SizeType size_type;
  typedef Kokkos::ArithTraits<typename YMV::non_const_value_type> ATS;

  XMV m_x;
  YMV m_y;
  const typename XMV::non_const_value_type m_a;
  const typename YMV::non_const_value_type m_b;

  Axpby_MV_Unroll_Functor(const XMV& X, const YMV& Y, const typename XMV::non_const_value_type& a,
                          const typename YMV::non_const_value_type& b, const SizeType /* startingColumn */)
      : m_x(X), m_y(Y), m_a(a), m_b(b) {
    static_assert(Kokkos::is_view<XMV>::value,
                  "KokkosBlas::Impl::Axpby_MV_Unroll_Functor(ABscalars)"
                  ": X is not a Kokkos::View.");
    static_assert(Kokkos::is_view<YMV>::value,
                  "KokkosBlas::Impl::Axpby_MV_Unroll_Functor(ABscalars)"
                  ": Y is not a Kokkos::View.");
    static_assert(std::is_same<typename YMV::value_type, typename YMV::non_const_value_type>::value,
                  "KokkosBlas::Impl::Axpby_MV_Unroll_Functor(ABscalars)"
                  ": Y must be nonconst, since it is an output argument"
                  " and we have to be able to write to its entries.");
    static_assert((int)YMV::rank == (int)XMV::rank,
                  "KokkosBlas::Impl::Axpby_MV_Unroll_Functor(ABscalars)"
                  ": X and Y must have the same rank.");
    static_assert(YMV::rank == 2,
                  "KokkosBlas::Impl::Axpby_MV_Unroll_Functor(ABscalars)"
                  ": XMV and YMV must have rank 2.");
    static_assert((-1 <= scalar_x) && (scalar_x <= 2) && (-1 <= scalar_y) && (scalar_y <= 2),
                  "KokkosBlas::Impl::Axpby_MV_Unroll_Functor(ABscalars)"
                  ": scalar_x and/or scalar_y are out of range.");
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_type& i) const {
    // scalar_x and scalar_y are compile-time constants (since they
    // are template parameters), so the compiler should evaluate these
    // branches at compile time.

    // **************************************************************
    // Possibilities with 'scalar_x == 0'
    // **************************************************************
    if constexpr (scalar_x == 0) {
      if constexpr (scalar_y == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(i, k) = ATS::zero();
        }
      } else if constexpr (scalar_y == -1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(i, k) = -m_y(i, k);
        }
      } else if constexpr (scalar_y == 1) {
        // Nothing to do: Y(i,j) := Y(i,j)
      } else if constexpr (scalar_y == 2) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(i, k) = m_b * m_y(i, k);
        }
      }
    }
    // **************************************************************
    // Possibilities with 'scalar_x == -1'
    // **************************************************************
    else if constexpr (scalar_x == -1) {
      if constexpr (scalar_y == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(i, k) = -m_x(i, k);
        }
      } else if constexpr (scalar_y == -1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(i, k) = -m_x(i, k) - m_y(i, k);
        }
      } else if constexpr (scalar_y == 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(i, k) = -m_x(i, k) + m_y(i, k);
        }
      } else if constexpr (scalar_y == 2) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(i, k) = -m_x(i, k) + m_b * m_y(i, k);
        }
      }
    }
    // **************************************************************
    // Possibilities with 'scalar_x == 1'
    // **************************************************************
    else if constexpr (scalar_x == 1) {
      if constexpr (scalar_y == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(i, k) = m_x(i, k);
        }
      } else if constexpr (scalar_y == -1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(i, k) = m_x(i, k) - m_y(i, k);
        }
      } else if constexpr (scalar_y == 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(i, k) = m_x(i, k) + m_y(i, k);
        }
      } else if constexpr (scalar_y == 2) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(i, k) = m_x(i, k) + m_b * m_y(i, k);
        }
      }
    }
    // **************************************************************
    // Possibilities with 'scalar_x == 2'
    // **************************************************************
    else if constexpr (scalar_x == 2) {
      if constexpr (scalar_y == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(i, k) = m_a * m_x(i, k);
        }
      } else if constexpr (scalar_y == -1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(i, k) = m_a * m_x(i, k) - m_y(i, k);
        }
      } else if constexpr (scalar_y == 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(i, k) = m_a * m_x(i, k) + m_y(i, k);
        }
      } else if constexpr (scalar_y == 2) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(i, k) = m_a * m_x(i, k) + m_b * m_y(i, k);
        }
      }
    }
  }
};

// Invoke the unrolled multivector functor that computes any of the
// following:
//
// 1. Y(i,j) = a*X(i,j) + b*Y(i,j) for a,b in -1,0,1
// 2. Y(i,j) = av(j)*X(i,j) + b*Y(i,j) for b in -1,0,1
// 3. Y(i,j) = a*X(i,j) + b*Y(i,j) for a in -1,0,1
// 4. Y(i,j) = av(j)*X(i,j) + bv(j)*Y(i,j)
//
// scalar_x and scalar_y come in as integers.  The values -1, 0, and 1
// correspond to the literal values of the coefficients.  The value 2
// tells the functor to use the corresponding vector of coefficients:
// - scalar_x == 2 means use av, otherwise ignore av;
// - scalar_y == 2 means use bv, otherwise ignore bv.
//
// Any literal coefficient of zero has BLAS semantics of ignoring the
// corresponding (multi)vector entry.  This does NOT apply to
// coefficients in av and bv vectors, if they are used.
//
// Either av and bv are both 1-D Views, or av and bv are both scalars.
template <class execution_space, class AV, class XMV, class BV, class YMV, int UNROLL, class SizeType>
void Axpby_MV_Unrolled(const execution_space& space, const AV& av, const XMV& x, const BV& bv, const YMV& y,
                       const SizeType startingColumn, int scalar_x = 2, int scalar_y = 2) {
  static_assert(Kokkos::is_view<XMV>::value,
                "KokkosBlas::Impl::Axpby_MV_Unrolled()"
                ": X is not a Kokkos::View.");
  static_assert(Kokkos::is_view<YMV>::value,
                "KokkosBlas::Impl::Axpby_MV_Unrolled()"
                ": Y is not a Kokkos::View.");
  static_assert(std::is_same<typename YMV::value_type, typename YMV::non_const_value_type>::value,
                "KokkosBlas::Impl::Axpby_MV_Unrolled()"
                ": Y must be nonconst, since it is an output argument"
                " and we have to be able to write to its entries.");
  static_assert((int)YMV::rank == (int)XMV::rank,
                "KokkosBlas::Impl::Axpby_MV_Unrolled()"
                ": X and Y must have the same rank.");
  static_assert(YMV::rank == 2,
                "KokkosBlas::Impl::Axpby_MV_Unrolled()"
                ": XMV and YMV must have rank 2.");
  if ((-1 <= scalar_x) && (scalar_x <= 2) && (-1 <= scalar_y) && (scalar_y <= 2)) {
    // Ok
  } else {
    KokkosKernels::Impl::throw_runtime_exception(
        "KokkosBlas::Impl::Axpby_MV_Unrolled()"
        ": scalar_x and/or scalar_y are out of range.");
  }

  const SizeType numRows = x.extent(0);
  Kokkos::RangePolicy<execution_space, SizeType> policy(space, 0, numRows);

  // ****************************************************************
  // Possibilities with 'scalar_x == 0'
  // ****************************************************************
  if (scalar_x == 0) {
    if (scalar_y == 0) {
      Axpby_MV_Unroll_Functor<AV, XMV, BV, YMV, 0, 0, UNROLL, SizeType> op(x, y, av, bv, startingColumn);
      Kokkos::parallel_for("KokkosBlas::Axpby::MV::S0", policy, op);
    } else if (scalar_y == -1) {
      Axpby_MV_Unroll_Functor<AV, XMV, BV, YMV, 0, -1, UNROLL, SizeType> op(x, y, av, bv, startingColumn);
      Kokkos::parallel_for("KokkosBlas::Axpby::MV::S1", policy, op);
    } else if (scalar_y == 1) {
      Axpby_MV_Unroll_Functor<AV, XMV, BV, YMV, 0, 1, UNROLL, SizeType> op(x, y, av, bv, startingColumn);
      Kokkos::parallel_for("KokkosBlas::Axpby::MV::S2", policy, op);
    } else if (scalar_y == 2) {
      Axpby_MV_Unroll_Functor<AV, XMV, BV, YMV, 0, 2, UNROLL, SizeType> op(x, y, av, bv, startingColumn);
      Kokkos::parallel_for("KokkosBlas::Axpby::MV::S3", policy, op);
    }
  }
  // ****************************************************************
  // Possibilities with 'scalar_x == -1'
  // ****************************************************************
  else if (scalar_x == -1) {
    if (scalar_y == 0) {
      Axpby_MV_Unroll_Functor<AV, XMV, BV, YMV, -1, 0, UNROLL, SizeType> op(x, y, av, bv, startingColumn);
      Kokkos::parallel_for("KokkosBlas::Axpby::MV::S4", policy, op);
    } else if (scalar_y == -1) {
      Axpby_MV_Unroll_Functor<AV, XMV, BV, YMV, -1, -1, UNROLL, SizeType> op(x, y, av, bv, startingColumn);
      Kokkos::parallel_for("KokkosBlas::Axpby::MV::S5", policy, op);
    } else if (scalar_y == 1) {
      Axpby_MV_Unroll_Functor<AV, XMV, BV, YMV, -1, 1, UNROLL, SizeType> op(x, y, av, bv, startingColumn);
      Kokkos::parallel_for("KokkosBlas::Axpby::MV::S6", policy, op);
    } else if (scalar_y == 2) {
      Axpby_MV_Unroll_Functor<AV, XMV, BV, YMV, -1, 2, UNROLL, SizeType> op(x, y, av, bv, startingColumn);
      Kokkos::parallel_for("KokkosBlas::Axpby::MV::S7", policy, op);
    }
  }
  // ****************************************************************
  // Possibilities with 'scalar_x == 1'
  // ****************************************************************
  else if (scalar_x == 1) {
    if (scalar_y == 0) {
      Axpby_MV_Unroll_Functor<AV, XMV, BV, YMV, 1, 0, UNROLL, SizeType> op(x, y, av, bv, startingColumn);
      Kokkos::parallel_for("KokkosBlas::Axpby::MV::S8", policy, op);
    } else if (scalar_y == -1) {
      Axpby_MV_Unroll_Functor<AV, XMV, BV, YMV, 1, -1, UNROLL, SizeType> op(x, y, av, bv, startingColumn);
      Kokkos::parallel_for("KokkosBlas::Axpby::MV::S9", policy, op);
    } else if (scalar_y == 1) {
      Axpby_MV_Unroll_Functor<AV, XMV, BV, YMV, 1, 1, UNROLL, SizeType> op(x, y, av, bv, startingColumn);
      Kokkos::parallel_for("KokkosBlas::Axpby::MV::S10", policy, op);
    } else if (scalar_y == 2) {
      Axpby_MV_Unroll_Functor<AV, XMV, BV, YMV, 1, 2, UNROLL, SizeType> op(x, y, av, bv, startingColumn);
      Kokkos::parallel_for("KokkosBlas::Axpby::MV::S11", policy, op);
    }
  }
  // ****************************************************************
  // Possibilities with 'scalar_x == 2'
  // ****************************************************************
  else if (scalar_x == 2) {
    if (scalar_y == 0) {
      Axpby_MV_Unroll_Functor<AV, XMV, BV, YMV, 2, 0, UNROLL, SizeType> op(x, y, av, bv, startingColumn);
      Kokkos::parallel_for("KokkosBlas::Axpby::MV::S12", policy, op);
    } else if (scalar_y == -1) {
      Axpby_MV_Unroll_Functor<AV, XMV, BV, YMV, 2, -1, UNROLL, SizeType> op(x, y, av, bv, startingColumn);
      Kokkos::parallel_for("KokkosBlas::Axpby::MV::S13", policy, op);
    } else if (scalar_y == 1) {
      Axpby_MV_Unroll_Functor<AV, XMV, BV, YMV, 2, 1, UNROLL, SizeType> op(x, y, av, bv, startingColumn);
      Kokkos::parallel_for("KokkosBlas::Axpby::MV::S14", policy, op);
    } else if (scalar_y == 2) {
      Axpby_MV_Unroll_Functor<AV, XMV, BV, YMV, 2, 2, UNROLL, SizeType> op(x, y, av, bv, startingColumn);
      Kokkos::parallel_for("KokkosBlas::Axpby::MV::S15", policy, op);
    }
  }
}

// Invoke the "generic" (not unrolled) multivector functor that
// computes any of the following:
//
// 1. Y(i,j) = a*X(i,j) + b*Y(i,j) for a,b in -1,0,1
// 2. Y(i,j) = av(j)*X(i,j) + b*Y(i,j) for b in -1,0,1
// 3. Y(i,j) = a*X(i,j) + b*Y(i,j) for a in -1,0,1
// 4. Y(i,j) = av(j)*X(i,j) + bv(j)*Y(i,j)
//
// scalar_x and scalar_y come in as integers.  The values -1, 0, and 1
// correspond to the literal values of the coefficients.  The value 2
// tells the functor to use the corresponding vector of coefficients:
// - scalar_x == 2 means use av, otherwise ignore av;
// - scalar_y == 2 means use bv, otherwise ignore bv.
//
// Any literal coefficient of zero has BLAS semantics of ignoring the
// corresponding (multi)vector entry.  This does NOT apply to
// coefficients in av and bv vectors, if they are used.
//
// Either av and bv are both 1-D Views, or av and bv are both scalars.
template <class execution_space, class AV, class XMV, class BV, class YMV, class SizeType>
void Axpby_MV_Generic(const execution_space& space, const AV& av, const XMV& x, const BV& bv, const YMV& y,
                      int scalar_x = 2, int scalar_y = 2) {
  static_assert(Kokkos::is_view<XMV>::value,
                "KokkosBlas::Impl::Axpby_MV_Generic()"
                ": X is not a Kokkos::View.");
  static_assert(Kokkos::is_view<YMV>::value,
                "KokkosBlas::Impl::Axpby_MV_Generic()"
                ": Y is not a Kokkos::View.");
  static_assert(std::is_same<typename YMV::value_type, typename YMV::non_const_value_type>::value,
                "KokkosBlas::Impl::Axpby_MV_Generic()"
                ": Y must be nonconst, since it is an output argument"
                " and we have to be able to write to its entries.");
  static_assert((int)YMV::rank == (int)XMV::rank,
                "KokkosBlas::Impl::Axpby_MV_Generic()"
                ": X and Y must have the same rank.");
  static_assert(YMV::rank == 2,
                "KokkosBlas::Impl::Axpby_MV_Generic()"
                ": XMV and YMV must have rank 2.");
  if ((-1 <= scalar_x) && (scalar_x <= 2) && (-1 <= scalar_y) && (scalar_y <= 2)) {
    // Ok
  } else {
    KokkosKernels::Impl::throw_runtime_exception(
        "KokkosBlas::Impl::Axpby_MV_Generic()"
        ": scalar_x and/or scalar_y are out of range.");
  }

  const SizeType numRows = x.extent(0);
  Kokkos::RangePolicy<execution_space, SizeType> policy(space, 0, numRows);

  // ****************************************************************
  // Possibilities with 'scalar_x == 0'
  // ****************************************************************
  if (scalar_x == 0) {
    if (scalar_y == 0) {
      Axpby_MV_Functor<AV, XMV, BV, YMV, 0, 0, SizeType> op(x, y, av, bv);
      Kokkos::parallel_for("KokkosBlas::Axpby::MV::S16", policy, op);
    } else if (scalar_y == -1) {
      Axpby_MV_Functor<AV, XMV, BV, YMV, 0, -1, SizeType> op(x, y, av, bv);
      Kokkos::parallel_for("KokkosBlas::Axpby::MV::S17", policy, op);
    } else if (scalar_y == 1) {
      Axpby_MV_Functor<AV, XMV, BV, YMV, 0, 1, SizeType> op(x, y, av, bv);
      Kokkos::parallel_for("KokkosBlas::Axpby::MV::S18", policy, op);
    } else if (scalar_y == 2) {
      Axpby_MV_Functor<AV, XMV, BV, YMV, 0, 2, SizeType> op(x, y, av, bv);
      Kokkos::parallel_for("KokkosBlas::Axpby::MV::S19", policy, op);
    }
  }
  // ****************************************************************
  // Possibilities with 'scalar_x == -1'
  // ****************************************************************
  else if (scalar_x == -1) {
    if (scalar_y == 0) {
      Axpby_MV_Functor<AV, XMV, BV, YMV, -1, 0, SizeType> op(x, y, av, bv);
      Kokkos::parallel_for("KokkosBlas::Axpby::MV::S20", policy, op);
    } else if (scalar_y == -1) {
      Axpby_MV_Functor<AV, XMV, BV, YMV, -1, -1, SizeType> op(x, y, av, bv);
      Kokkos::parallel_for("KokkosBlas::Axpby::MV::S21", policy, op);
    } else if (scalar_y == 1) {
      Axpby_MV_Functor<AV, XMV, BV, YMV, -1, 1, SizeType> op(x, y, av, bv);
      Kokkos::parallel_for("KokkosBlas::Axpby::MV::S22", policy, op);
    } else if (scalar_y == 2) {
      Axpby_MV_Functor<AV, XMV, BV, YMV, -1, 2, SizeType> op(x, y, av, bv);
      Kokkos::parallel_for("KokkosBlas::Axpby::MV::S23", policy, op);
    }
  }
  // ****************************************************************
  // Possibilities with 'scalar_x == 1'
  // ****************************************************************
  else if (scalar_x == 1) {
    if (scalar_y == 0) {
      Axpby_MV_Functor<AV, XMV, BV, YMV, 1, 0, SizeType> op(x, y, av, bv);
      Kokkos::parallel_for("KokkosBlas::Axpby::MV::S24", policy, op);
    } else if (scalar_y == -1) {
      Axpby_MV_Functor<AV, XMV, BV, YMV, 1, -1, SizeType> op(x, y, av, bv);
      Kokkos::parallel_for("KokkosBlas::Axpby::MV::S25", policy, op);
    } else if (scalar_y == 1) {
      Axpby_MV_Functor<AV, XMV, BV, YMV, 1, 1, SizeType> op(x, y, av, bv);
      Kokkos::parallel_for("KokkosBlas::Axpby::MV::S26", policy, op);
    } else if (scalar_y == 2) {
      Axpby_MV_Functor<AV, XMV, BV, YMV, 1, 2, SizeType> op(x, y, av, bv);
      Kokkos::parallel_for("KokkosBlas::Axpby::MV::S27", policy, op);
    }
  }
  // ****************************************************************
  // Possibilities with 'scalar_x == 2'
  // ****************************************************************
  else if (scalar_x == 2) {
    if (scalar_y == 0) {
      Axpby_MV_Functor<AV, XMV, BV, YMV, 2, 0, SizeType> op(x, y, av, bv);
      Kokkos::parallel_for("KokkosBlas::Axpby::MV::S28", policy, op);
    } else if (scalar_y == -1) {
      Axpby_MV_Functor<AV, XMV, BV, YMV, 2, -1, SizeType> op(x, y, av, bv);
      Kokkos::parallel_for("KokkosBlas::Axpby::MV::S29", policy, op);
    } else if (scalar_y == 1) {
      Axpby_MV_Functor<AV, XMV, BV, YMV, 2, 1, SizeType> op(x, y, av, bv);
      Kokkos::parallel_for("KokkosBlas::Axpby::MV::S30", policy, op);
    } else if (scalar_y == 2) {
      Axpby_MV_Functor<AV, XMV, BV, YMV, 2, 2, SizeType> op(x, y, av, bv);
      Kokkos::parallel_for("KokkosBlas::Axpby::MV::S31", policy, op);
    }
  }
}

// Compute any of the following, in a way optimized for X and Y
// being LayoutLeft:
//
// 1. Y(i,j) = a*X(i,j) + b*Y(i,j) for a,b in -1,0,1
// 2. Y(i,j) = av(j)*X(i,j) + b*Y(i,j) for b in -1,0,1
// 3. Y(i,j) = a*X(i,j) + b*Y(i,j) for a in -1,0,1
// 4. Y(i,j) = av(j)*X(i,j) + bv(j)*Y(i,j)
//
// scalar_x and scalar_y come in as integers.  The values -1, 0, and 1
// correspond to the literal values of the coefficients.  The value 2
// tells the functor to use the corresponding vector of coefficients:
// - scalar_x == 2 means use av, otherwise ignore av;
// - scalar_y == 2 means use bv, otherwise ignore bv.
//
// Any literal coefficient of zero has BLAS semantics of ignoring the
// corresponding (multi)vector entry.  This does NOT apply to
// coefficients in av and bv vectors, if they are used.
//
// Either av and bv are both 1-D Views, or av and bv are both scalars.
template <class execution_space, class AV, class XMV, class BV, class YMV, class SizeType>
struct Axpby_MV_Invoke_Left {
  static void run(const execution_space& space, const AV& av, const XMV& x, const BV& bv, const YMV& y,
                  int scalar_x = 2, int scalar_y = 2) {
    static_assert(Kokkos::is_view<XMV>::value,
                  "KokkosBlas::Impl::Axpby_MV_Invoke_Left::run()"
                  ": X is not a Kokkos::View.");
    static_assert(Kokkos::is_view<YMV>::value,
                  "KokkosBlas::Impl::Axpby_MV_Invoke_Left::run()"
                  ": Y is not a Kokkos::View.");
    static_assert(std::is_same<typename YMV::value_type, typename YMV::non_const_value_type>::value,
                  "KokkosBlas::Impl::Axpby_MV_Invoke_Left::run()"
                  ": Y must be nonconst, since it is an output argument"
                  " and we have to be able to write to its entries.");
    static_assert((int)YMV::rank == (int)XMV::rank,
                  "KokkosBlas::Impl::Axpby_MV_Invoke_Left::run()"
                  ": X and Y must have the same rank.");
    static_assert(YMV::rank == 2,
                  "KokkosBlas::Impl::Axpby_MV_Invoke_Left::run()"
                  ": X and Y must have rank 2.");
    if ((-1 <= scalar_x) && (scalar_x <= 2) && (-1 <= scalar_y) && (scalar_y <= 2)) {
      // Ok
    } else {
      KokkosKernels::Impl::throw_runtime_exception(
          "KokkosBlas::Impl::Axpby_MV_Invoke_Left::run()"
          ": scalar_x and/or scalar_y are out of range.");
    }

    const SizeType numCols = x.extent(1);

    // Strip-mine by 8, then 4.  After that, do one column at a time.
    // We limit the number of strip-mine values in order to keep down
    // the amount of code to compile.
    SizeType j = 0;
    for (; j + 8 <= numCols; j += 8) {
      XMV X_cur = Kokkos::subview(x, Kokkos::ALL(), std::make_pair(j, j + 8));
      YMV Y_cur = Kokkos::subview(y, Kokkos::ALL(), std::make_pair(j, j + 8));

      // Passing in the starting column index lets the functor take
      // subviews of av and bv, if they are Views.  If they are scalars,
      // the functor doesn't have to do anything to them.
      Axpby_MV_Unrolled<execution_space, AV, XMV, BV, YMV, 8, SizeType>(space, av, X_cur, bv, Y_cur, j, scalar_x,
                                                                        scalar_y);
    }
    for (; j + 4 <= numCols; j += 4) {
      XMV X_cur = Kokkos::subview(x, Kokkos::ALL(), std::make_pair(j, j + 4));
      YMV Y_cur = Kokkos::subview(y, Kokkos::ALL(), std::make_pair(j, j + 4));

      // Passing in the starting column index lets the functor take
      // subviews of av and bv, if they are Views.  If they are scalars,
      // the functor doesn't have to do anything to them.
      Axpby_MV_Unrolled<execution_space, AV, XMV, BV, YMV, 4, SizeType>(space, av, X_cur, bv, Y_cur, j, scalar_x,
                                                                        scalar_y);
    }
    for (; j < numCols; ++j) {
      auto x_cur = Kokkos::subview(x, Kokkos::ALL(), j);
      auto y_cur = Kokkos::subview(y, Kokkos::ALL(), j);

      // Passing in the starting column index lets the functor take
      // subviews of av and bv, if they are Views.  If they are scalars,
      // the functor doesn't have to do anything to them.
      typedef decltype(x_cur) XV;
      typedef decltype(y_cur) YV;
      Axpby_Generic<execution_space, AV, XV, BV, YV, SizeType>(space, av, x_cur, bv, y_cur, j, scalar_x, scalar_y);
    }
  }
};

// Compute any of the following, in a way optimized for X and Y
// being LayoutRight:
//
// 1. Y(i,j) = a*X(i,j) + b*Y(i,j) for a,b in -1,0,1
// 2. Y(i,j) = av(j)*X(i,j) + b*Y(i,j) for b in -1,0,1
// 3. Y(i,j) = a*X(i,j) + bv(j)*Y(i,j) for a in -1,0,1
// 4. Y(i,j) = av(j)*X(i,j) + bv(j)*Y(i,j)
//
// scalar_x and scalar_y come in as integers.  The values -1, 0, and 1
// correspond to the literal values of the coefficients.  The value 2
// tells the functor to use the corresponding vector of coefficients:
// - scalar_x == 2 means use av, otherwise ignore av;
// - scalar_y == 2 means use bv, otherwise ignore bv.
//
// Any literal coefficient of zero has BLAS semantics of ignoring the
// corresponding (multi)vector entry.  This does NOT apply to
// coefficients in av and bv vectors, if they are used.
//
// Either av and bv are both 1-D Views, or av and bv are both scalars.
template <class execution_space, class AV, class XMV, class BV, class YMV, class SizeType>
struct Axpby_MV_Invoke_Right {
  static void run(const execution_space& space, const AV& av, const XMV& x, const BV& bv, const YMV& y,
                  int scalar_x = 2, int scalar_y = 2) {
    static_assert(Kokkos::is_view<XMV>::value,
                  "KokkosBlas::Impl::Axpby_MV_Invoke_Right::run()"
                  ": X is not a Kokkos::View.");
    static_assert(Kokkos::is_view<YMV>::value,
                  "KokkosBlas::Impl::Axpby_MV_Invoke_Right::run()"
                  ": Y is not a Kokkos::View.");
    static_assert(std::is_same<typename YMV::value_type, typename YMV::non_const_value_type>::value,
                  "KokkosBlas::Impl::Axpby_MV_Invoke_Right::run()"
                  ": Y must be nonconst, since it is an output argument"
                  " and we have to be able to write to its entries.");
    static_assert((int)YMV::rank == (int)XMV::rank,
                  "KokkosBlas::Impl::Axpby_MV_Invoke_Right::run()"
                  ": X and Y must have the same rank.");
    static_assert(YMV::rank == 2,
                  "KokkosBlas::Impl::Axpby_MV_Invoke_Right::run()"
                  ": X and Y must have rank 2.");
    if ((-1 <= scalar_x) && (scalar_x <= 2) && (-1 <= scalar_y) && (scalar_y <= 2)) {
      // Ok
    } else {
      KokkosKernels::Impl::throw_runtime_exception(
          "KokkosBlas::Impl::Axpby_MV_Invoke_Right::run()"
          ": scalar_x and/or scalar_y are out of range.");
    }

    const SizeType numCols = x.extent(1);
    if (numCols == 1) {
      auto x_0 = Kokkos::subview(x, Kokkos::ALL(), 0);
      auto y_0 = Kokkos::subview(y, Kokkos::ALL(), 0);
      typedef decltype(x_0) XV;
      typedef decltype(y_0) YV;
      Axpby_Generic<execution_space, AV, XV, BV, YV, SizeType>(space, av, x_0, bv, y_0, 0, scalar_x, scalar_y);
    } else {
      Axpby_MV_Generic<execution_space, AV, XMV, BV, YMV, SizeType>(space, av, x, bv, y, scalar_x, scalar_y);
    }
  }
};

}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSBLAS1_AXPBY_MV_IMPL_HPP_
