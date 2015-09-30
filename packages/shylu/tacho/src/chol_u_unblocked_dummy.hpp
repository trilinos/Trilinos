#pragma once
#ifndef __CHOL_U_UNBLOCKED_DUMMY_HPP__
#define __CHOL_U_UNBLOCKED_DUMMY_HPP__

/// \file chol_u_unblocked_dummy.hpp
/// \brief Test code for data parallel interface
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "util.hpp"
#include "partition.hpp"

namespace Tacho {

  using namespace std;

  template<template<int> class ControlType,
           typename ParallelForType,
           typename CrsExecViewType>
  class FunctionChol<Uplo::Upper,AlgoChol::Dummy,
                     ControlType,
                     ParallelForType,
                     CrsExecViewType> {
  private:
    int _r_val;

  public:
    KOKKOS_INLINE_FUNCTION
    int getReturnValue() const {
      return _r_val;
    }

    KOKKOS_INLINE_FUNCTION
    FunctionChol(typename CrsExecViewType::policy_type &policy,
                 const typename CrsExecViewType::policy_type::member_type &member,
                 CrsExecViewType &A) {
      typedef typename CrsExecViewType::value_type        value_type;
      typedef typename CrsExecViewType::ordinal_type      ordinal_type;
      typedef typename CrsExecViewType::row_view_type     row_view_type;
      typedef typename CrsExecViewType::team_factory_type team_factory_type;
      
      // row_view_type r1t, r2t;
      
      for (ordinal_type k=0;k<A.NumRows();++k) {
        //r1t.setView(A, k);
        row_view_type &r1t = A.RowView(k);
        
        // extract diagonal from alpha11
        value_type &alpha = r1t.Value(0);
        const ordinal_type nnz_r1t = r1t.NumNonZeros();
        
        if (nnz_r1t) {
          // inverse scale
          ParallelForType(team_factory_type::createThreadLoopRegion(member, 1, nnz_r1t),
                          [&](const ordinal_type j) {
                            r1t.Value(j) /= alpha;
                          });
        }
      }
      _r_val = 0;
    }
  };
}

#endif
