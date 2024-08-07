// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Zoltan2_AlgMetis.hpp
    \brief The ND ordering algorithm uses Metis.
*/

#ifndef _ZOLTAN2_ALGMETIS_HPP_
#define _ZOLTAN2_ALGMETIS_HPP_

#include <Zoltan2_Algorithm.hpp>
#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_OrderingSolution.hpp>
#include <Zoltan2_TPLTraits.hpp>
#ifdef HAVE_ZOLTAN2_METIS
#include "metis.h"
#endif

namespace Zoltan2{

template <typename Adapter>
class AlgMetis : public Algorithm<Adapter>
{
    private:

    const RCP<const typename Adapter::base_adapter_t> adapter;
    const RCP<Teuchos::ParameterList> pl;
    const RCP<const Teuchos::Comm<int> > comm;
    RCP<const Environment> env;
    modelFlag_t graphFlags;

    public:

    AlgMetis(
      const RCP<const typename Adapter::base_adapter_t> &adapter__,
      const RCP<Teuchos::ParameterList> &pl__,
      const RCP<const Teuchos::Comm<int> > &comm__,
      RCP<const Environment> &env__,
      const modelFlag_t &graphFlags__
    ) : adapter(adapter__), pl(pl__), comm(comm__), env(env__), graphFlags(graphFlags__)
    { }

    int globalOrder(
      const RCP<GlobalOrderingSolution<typename Adapter::gno_t> > &/* solution */) {
        throw std::logic_error("AlgMetis does not yet support global ordering.");
    }

    int localOrder(
      const RCP<LocalOrderingSolution<typename Adapter::lno_t> > &solution)
    {
#ifndef HAVE_ZOLTAN2_METIS
      (void)solution; // remove unused parameter warning
  throw std::runtime_error(
        "BUILD ERROR: Metis requested but not compiled into Zoltan2.\n"
        "Please set CMake flag Zoltan2_ENABLE_METIS:BOOL=ON.");
#else
      typedef typename Adapter::gno_t gno_t;
      typedef typename Adapter::lno_t lno_t;
      typedef typename Adapter::offset_t offset_t;
      typedef typename Adapter::scalar_t scalar_t;

      int ierr= 0;

      // Get EdgeList
      const auto model = rcp(new GraphModel<Adapter>(adapter, env, comm, graphFlags));
      const size_t nVtx = model->getLocalNumVertices();
      const size_t nNnz = model->getLocalNumEdges();
      lno_t *perm = (lno_t *) (solution->getPermutationRCP().getRawPtr());

      if (nVtx > 0 && nNnz > 0) {
        ArrayView<const gno_t> edgeIds;
        ArrayView<const offset_t> offsets;
        ArrayView<StridedData<lno_t, scalar_t> > wgts; // wgts are ignored in NodeND
        model->getEdgeList(edgeIds, offsets, wgts);

        // Prepare for calling metis
        using Zoltan2OffsetView = typename Kokkos::View<offset_t*, Kokkos::HostSpace>;
        using Zoltan2EdgeView = typename Kokkos::View<gno_t*, Kokkos::HostSpace>;
        Zoltan2OffsetView zoltan2_rowptr (const_cast<offset_t*>(offsets.data()), nVtx+1);
        Zoltan2EdgeView   zoltan2_colidx (const_cast<gno_t*>(edgeIds.data()), nNnz);

        using MetisIdxView = typename Kokkos::View<idx_t*, Kokkos::HostSpace>;
        MetisIdxView metis_rowptr;
        MetisIdxView metis_colidx;

        // Symmetrize (always for now)
        KokkosKernels::Impl::symmetrize_graph_symbolic_hashmap<
          Zoltan2OffsetView, Zoltan2EdgeView, MetisIdxView, MetisIdxView, Kokkos::HostSpace::execution_space>
          (nVtx, zoltan2_rowptr, zoltan2_colidx, metis_rowptr, metis_colidx);

        // Remove diagonals
        idx_t metis_nVtx=0;
        TPL_Traits<idx_t, size_t>::ASSIGN(metis_nVtx, nVtx);

        idx_t nnz = metis_rowptr(0);
        idx_t old_nnz = nnz;
        for (idx_t i = 0; i < metis_nVtx; i++) {
          for (idx_t k = old_nnz; k < metis_rowptr(i+1); k++) {
            if (metis_colidx(k) != i) {
              metis_colidx(nnz) = metis_colidx(k);
              nnz++;
            }
          }
          old_nnz = metis_rowptr(i+1);
          metis_rowptr(i+1) = nnz;
        }

        // Allocate Metis perm/iperm
        idx_t *metis_perm = new idx_t[nVtx];
        idx_t *metis_iperm = new idx_t[nVtx];

        // Call metis
        int info = METIS_NodeND(&metis_nVtx, metis_rowptr.data(), metis_colidx.data(),
                                NULL, NULL, metis_perm, metis_iperm);
        if (METIS_OK != info) {
          throw std::runtime_error(std::string("METIS_NodeND returned info = " + info));
        }

        // Copy result back
        for (size_t i = 0; i < nVtx; i++)
          TPL_Traits<lno_t, idx_t>::ASSIGN(perm[i], metis_perm[i]);

        delete [] metis_iperm;
        delete [] metis_perm;
      } else {
        for (size_t i = 0; i < nVtx; i++)
          perm[i] = i;
      }

      solution->setHavePerm(true);
      return ierr;
#endif
    }
};

}



#endif
