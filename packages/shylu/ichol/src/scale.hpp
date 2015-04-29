#pragma once
#ifndef __SCALE_HPP__
#define __SCALE_HPP__

/// \file scale.hpp
/// \brief Scaling sparse matrix.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Example {

  using namespace std;

  template<typename T> struct ScaleTraits {
    typedef T scale_type;
    static T one;
    static T zero;
  };

  // assume built-in types have appropriate type conversion
  template<typename T> T ScaleTraits<T>::one  = 1;
  template<typename T> T ScaleTraits<T>::zero = 0;

  template<typename ParallelForType,
           typename ScalarType,
           typename CrsExecViewType>
  KOKKOS_INLINE_FUNCTION
  int
  scaleCrsMatrix(const typename CrsExecViewType::policy_type::member_type &member,
                 const ScalarType &alpha, 
                 const CrsExecViewType &A) {
    typedef typename CrsExecViewType::ordinal_type  ordinal_type;
    typedef typename CrsExecViewType::value_type    value_type;
    typedef typename CrsExecViewType::row_view_type row_view_type;
    typedef typename CrsExecViewType::team_factory_type team_factory_type;

    if (alpha == ScaleTraits<value_type>::one) {
      // do nothing
    } else {
      ParallelForType(team_factory_type::createThreadLoopRegion(member, 0, A.NumRows()),
                      [&](const ordinal_type i) {
                        row_view_type row(A, i);
                        for (ordinal_type j=0;j<row.NumNonZeros();++j)
                          row.Value(j) *= alpha;
                      });
    }

    return 0;
  }

  template<typename ParallelForType,
           typename ScalarType,
           typename DenseExecViewType>
  KOKKOS_INLINE_FUNCTION
  int
  scaleDenseMatrix(const typename DenseExecViewType::policy_type::member_type &member,
                   const ScalarType &alpha, 
                   const DenseExecViewType &A) {
    typedef typename DenseExecViewType::ordinal_type  ordinal_type;
    typedef typename DenseExecViewType::value_type    value_type;
    typedef typename DenseExecViewType::team_factory_type team_factory_type;

    if (alpha == ScaleTraits<value_type>::one) {
      // do nothing
    } else {
      if (A.Base->ColStride() > A.Base->RowStride()) 
        ParallelForType(team_factory_type::createThreadLoopRegion(member, 0, A.NumCols()),
                        [&](const ordinal_type j) {
                          for (ordinal_type i=0;i<A.NumRows();++i)
                            A.Value(i, j) *= alpha;
                        });
      else
        ParallelForType(team_factory_type::createThreadLoopRegion(member, 0, A.NumRows()),
                        [&](const ordinal_type i) {
                          for (ordinal_type j=0;j<A.NumCols();++j)
                            A.Value(i, j) *= alpha;
                        });
    }

    return 0;
  }

}

#endif

