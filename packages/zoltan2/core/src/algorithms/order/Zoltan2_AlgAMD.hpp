// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Zoltan2_AlgAMD.hpp
    \brief The AMD ordering algorithm uses SuiteSparse.
*/

#ifndef _ZOLTAN2_ALGAMD_HPP_
#define _ZOLTAN2_ALGAMD_HPP_

#include <Zoltan2_Algorithm.hpp>
#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_OrderingSolution.hpp>
#include <Zoltan2_TPLTraits.hpp>


#ifdef HAVE_ZOLTAN2_AMD
#include "amd.h"
#ifdef SUITESPARSE_MAIN_VERSION
#if SUITESPARSE_MAIN_VERSION < 4
typedef UF_long SuiteSparse_long;
#endif
#endif
#endif



namespace Zoltan2{


#ifdef HAVE_ZOLTAN2_AMD
template <typename Ordinal>
class AMDTraits
{
    public:
    Ordinal order(Ordinal n, const Ordinal *Ap, const Ordinal *Ai,
                Ordinal *perm, double *control, double *info);
};

template <>
class AMDTraits<int>
{
    public:
    int order(int n, const int *Ap, const int *Ai, int *perm,
                double *control, double *info)
    {
        return (amd_order(n, Ap, Ai, perm, control, info));
    }
};

template <>
class AMDTraits<SuiteSparse_long>
{
    public:
    long order(SuiteSparse_long n, const SuiteSparse_long *Ap,
                const SuiteSparse_long *Ai, SuiteSparse_long *perm,
                double *control, double *info)
    {
        return (amd_l_order(n, Ap, Ai, perm, control, info));
    }
};

#endif

}


namespace Zoltan2{

template <typename Adapter>
class AlgAMD : public Algorithm<Adapter>
{
    private:

    const RCP<const typename Adapter::base_adapter_t> adapter;
    const RCP<Teuchos::ParameterList> pl;
    const RCP<const Teuchos::Comm<int> > comm;
    RCP<const Environment> env;
    modelFlag_t graphFlags;

    public:

    AlgAMD(
      const RCP<const typename Adapter::base_adapter_t> &adapter__,
      const RCP<Teuchos::ParameterList> &pl__,
      const RCP<const Teuchos::Comm<int> > &comm__,
      RCP<const Environment> &env__,
      const modelFlag_t &graphFlags__
    ) : adapter(adapter__), pl(pl__), comm(comm__), env(env__), graphFlags(graphFlags__)
    { }

    int globalOrder(
      const RCP<GlobalOrderingSolution<typename Adapter::gno_t> > &/* solution */) {
        throw std::logic_error("AlgAMD does not yet support global ordering.");
    }

    int localOrder(
      const RCP<LocalOrderingSolution<typename Adapter::lno_t> > &solution)
    {
#ifndef HAVE_ZOLTAN2_AMD
      (void)solution; // remove unused parameter warning
  throw std::runtime_error(
        "BUILD ERROR: AMD requested but not compiled into Zoltan2.\n"
        "Please set CMake flag Zoltan2_ENABLE_AMD:BOOL=ON.");
#else
      typedef typename Adapter::gno_t gno_t;
      typedef typename Adapter::lno_t lno_t;
      typedef typename Adapter::offset_t offset_t;
      typedef typename Adapter::scalar_t scalar_t;

      int ierr= 0;

      const auto model = rcp(new GraphModel<Adapter>(adapter, env, comm, graphFlags));
      const size_t nVtx = model->getLocalNumVertices();

      //cout << "Local num vertices" << nVtx << endl;
      ArrayView<const gno_t> edgeIds;
      ArrayView<const offset_t> offsets;
      ArrayView<StridedData<lno_t, scalar_t> > wgts;

      // wgts are ignored in AMD
      model->getEdgeList(edgeIds, offsets, wgts);

      // We will always call AMD with SuiteSparse_long
      AMDTraits<SuiteSparse_long> AMDobj;
      double Control[AMD_CONTROL];
      double Info[AMD_INFO];

      amd_defaults(Control);
      amd_control(Control);

      // We will use the lno_t for local ordering
      lno_t *perm;
      perm = (lno_t *) (solution->getPermutationRCP().getRawPtr());

      SuiteSparse_long *amd_offsets, *amd_edgeIds;
      TPL_Traits<SuiteSparse_long, const offset_t>::ASSIGN_ARRAY(&amd_offsets,
             offsets);
      // This might throw depending on how SuiteSparse was compiled
      // with long or long long and the size of both of them
      TPL_Traits<SuiteSparse_long, const gno_t>::ASSIGN_ARRAY(&amd_edgeIds,
             edgeIds);

      SuiteSparse_long amd_nVtx=0;
      TPL_Traits<SuiteSparse_long, size_t>::ASSIGN(amd_nVtx, nVtx);

      // Allocate a SuiteSparse_long perm
      SuiteSparse_long *amd_perm = new SuiteSparse_long[amd_nVtx];

      lno_t result = AMDobj.order(amd_nVtx, amd_offsets,
                             amd_edgeIds, amd_perm, Control, Info);

      if (result != AMD_OK && result != AMD_OK_BUT_JUMBLED)
          ierr = -1;

      // SR: This conversion might throw as we are going from SuiteSparse_long
      // to lno_t. Another option is to change local ordering solution
      // to use gno_t everywhere
      for (size_t i = 0; i < nVtx; i++)
          TPL_Traits<lno_t, SuiteSparse_long>::ASSIGN(perm[i], amd_perm[i]);

      // Delete copies
      TPL_Traits<SuiteSparse_long, const lno_t>::DELETE_ARRAY(&amd_offsets);
      TPL_Traits<SuiteSparse_long, const gno_t>::DELETE_ARRAY(&amd_edgeIds);
      delete [] amd_perm;

      solution->setHavePerm(true);
      return ierr;
#endif
    }
};

}



#endif
