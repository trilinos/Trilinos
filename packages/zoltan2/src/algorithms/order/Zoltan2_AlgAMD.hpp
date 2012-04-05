// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_AlgAMD.hpp
    \brief The AMD ordering algorithm.

  \todo  If Zoltan2 is not using int, int, double, we exit with an error.
           Do we want to copy the data to the right format for AMD?
*/

#ifndef _ZOLTAN2_ALGAMD_HPP_
#define _ZOLTAN2_ALGAMD_HPP_

#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_OrderingSolution.hpp>

#ifdef HAVE_ZOLTAN2_AMD
#include "amd.h"
#endif

template <typename Ordinal>
class AMDTraits
{
    public:
    Ordinal order(Ordinal n, const Ordinal *Ap, const Ordinal *Ai,
                Ordinal *perm, double *control, double *info);
};

#ifdef HAVE_ZOLTAN2_AMD
#ifdef HAVE_ZOLTAN2_INST_DOUBLE_INT_INT
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
class AMDTraits<long>
{
    public:
    long order(long n, const long *Ap, const long *Ai, long *perm,
                double *control, double *info)
    {
        return (amd_l_order(n, Ap, Ai, perm, control, info));
    }
};
#endif
#endif

namespace Zoltan2{

template <typename Adapter>
int AlgAMD(
  const RCP<GraphModel<Adapter> > &model,
  const RCP<OrderingSolution<typename Adapter::gid_t,
                             typename Adapter::lno_t> > &solution,
  const RCP<Teuchos::ParameterList> &pl,
  const RCP<const Teuchos::Comm<int> > &comm
) 
{
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::gid_t gid_t;
  typedef typename Adapter::scalar_t scalar_t;

  int ierr= 0;

  const size_t nVtx = model->getLocalNumVertices();

  cout << "Local num vertices" << nVtx << endl;
  ArrayView<const gno_t> edgeIds;
  ArrayView<const int> procIds;
  ArrayView<const lno_t> offsets;
  ArrayView<const scalar_t> wgts;

  //const size_t nEdgs = model->getEdgeList( edgeIds,
  //                      procIds, offsets, wgts);
  // TODO: Should use the local graph
  model->getEdgeList( edgeIds, procIds, offsets, wgts);


#ifdef HAVE_ZOLTAN2_AMD
#ifdef HAVE_ZOLTAN2_INST_DOUBLE_INT_INT
  cout << "AMD is enabled" << endl;
  AMDTraits<lno_t> AMDobj;
  double Control[AMD_CONTROL];
  double Info[AMD_INFO];

  amd_defaults(Control);
  amd_control(Control);

  lno_t *perm;
  perm = (lno_t *) (solution->getPermutationRCP().getRawPtr());

  lno_t result = AMDobj.order(nVtx, offsets.getRawPtr(),
                         edgeIds.getRawPtr(), perm, Control, Info);

  if (result != AMD_OK && result != AMD_OK_BUT_JUMBLED)
      ierr = -1; // TODO: Change return value to lno_t
#else
  cout << "Zoltan2 is not built with AMD's data types (int, int, double)." << endl;
  ierr = -2; // TODO: Need better error numbers ?
#endif

#else
  cout << "Zoltan2 is not built with AMD." << endl;
  ierr = -2; // TODO: Need better error numbers ?
#endif

  return ierr;
}

}

#endif
