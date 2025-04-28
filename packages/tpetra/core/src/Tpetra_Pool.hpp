// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// clang-format off

#ifndef TPETRA_POOL_HPP
#define TPETRA_POOL_HPP

#include "Kokkos_Core.hpp"
#include "Tpetra_MultiVector_decl.hpp"

namespace Tpetra
{

namespace impl
{
template<class DualViewType>
class DualViewPool
{
public:
  using dv_t = DualViewType;

  DualViewPool() {
    Kokkos::push_finalize_hook([this]() { this->availableDVs.clear(); });
  }

  template <typename... Sizes>
  Teuchos::RCP<dv_t> getDV(const Sizes & ...sizes){
    size_t total_size = dv_t::t_dev::required_allocation_size(sizes...);

    // Use lower_bound so that we can re-use a slightly larger allocation if it is available
    auto available_it = availableDVs.lower_bound(total_size);
    while(available_it != availableDVs.end() && available_it->second.empty()) {
      ++available_it;
    }
    if(available_it != availableDVs.end()) {
      auto & available = available_it->second;
      auto full_size_dv = available.back();
      available.pop_back();

      typename dv_t::t_dev mv_dev(full_size_dv.view_device().data(), sizes...);
      typename dv_t::t_host mv_host(full_size_dv.view_host().data(), sizes...);

      return Teuchos::rcpWithDealloc(new dv_t(mv_dev, mv_host), RCPDeleter{available, full_size_dv});
    }

    // No sufficiently large allocations were found so we need to create a new one.
    // Also remove the largest currently available allocation if there is one because it would be able
    // to use the allocation we are adding instead.
    auto available_rit = availableDVs.rbegin();
    while(available_rit != availableDVs.rend() && available_rit->second.empty()) {
      ++available_rit;
    }
    if(available_rit != availableDVs.rend()) {
      available_rit->second.pop_back();
    }
    dv_t dv("Tpetra::DualViewPool", sizes...);
    auto & available = availableDVs[total_size];
    return Teuchos::rcpWithDealloc(new dv_t(dv), RCPDeleter{available, dv});
  }

private:
  struct RCPDeleter
  {
    void free(dv_t * dv_ptr) {
      if(dv_ptr) {
        dv_pool.push_back(dv);
        delete dv_ptr;
      }
    }
    std::vector<dv_t> & dv_pool;
    dv_t dv;
  };
  std::map<size_t, std::vector<dv_t>> availableDVs;
};
}

template<class DualViewType>
inline impl::DualViewPool<DualViewType> & getPool()
{
  static impl::DualViewPool<DualViewType> static_pool;
  return static_pool;
}

template<class DualViewType, class... Sizes>
inline Teuchos::RCP<DualViewType> getDualViewFromPool(const Sizes & ...sizes)
{
  return getPool<DualViewType>().getDV(sizes...);
}

namespace impl
{
template<class Scalar, class LO, class GO, class Node>
class MultiVecPool
{
public:
  MultiVecPool() = default;

  using MV = ::Tpetra::MultiVector<Scalar, LO, GO, Node>;

  Teuchos::RCP<MV> getMV(const Teuchos::RCP<const ::Tpetra::Map<LO, GO, Node>> & map, const int numVecs) {
    const auto num_local_elems = map->getLocalNumElements();
    auto dv_rcp = getDualViewFromPool<dv_t>(num_local_elems, numVecs);

    return Teuchos::rcpWithDealloc(new MV(map, *dv_rcp), RCPDeleter{dv_rcp});
  }

private:
  using dv_t = typename MV::dual_view_type;
  struct RCPDeleter
  {
    void free(MV * mv_ptr) {
      if(mv_ptr) {
        delete mv_ptr;
        dv_rcp.reset();
      }
    }
    Teuchos::RCP<dv_t> dv_rcp;
  };
  std::map<size_t, std::vector<dv_t>> availableDVs;
};
}

template<class Scalar, class LO, class GO, class Node>
inline impl::MultiVecPool<Scalar, LO, GO, Node> & getPool()
{
  static impl::MultiVecPool<Scalar, LO, GO, Node> static_pool;
  return static_pool;
}

template<class Scalar, class LO, class GO, class Node>
inline Teuchos::RCP<::Tpetra::MultiVector<Scalar, LO, GO, Node>> getMultiVectorFromPool(const Teuchos::RCP<const ::Tpetra::Map<LO, GO, Node>> & map, const int numVecs)
{
  return getPool<Scalar, LO, GO, Node>().getMV(map, numVecs);
}

}

#endif