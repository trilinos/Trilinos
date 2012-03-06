// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_Metric.hpp
 *  \brief namespace methods to compute quality metrics
 *  \todo We need graph and hypergraph metrics.
 */

#ifndef ZOLTAN2_METRIC_HPP
#define ZOLTAN2_METRIC_HPP

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_Environment.hpp>

namespace Zoltan2{

/*! \brief Compute the global imbalance, one per weight dimension.
 *
 *  All processes in the communicator must call this.
 *
 *  \todo as written now this function uses buffers on the
 *    order of numGlobalParts.  This should be fixed eventually,
 *    when numGlobalParts may be very large.
 *
 *  \param  env    library environment information
 *
 *  \param  comm   the communicator
 *
 *  \param  numGlobalParts   the desired global number of parts
 *
 *  \param  partSizes has one array for each weight dimension.  If
 *     partSizes[w].size() is zero, then we assume uniform part sizes
 *     are desired for that weight dimension.  Otherwise partSizes[w] has
 *     size numGlobalParts, partSizes[w][p] is the part size desired for 
 *     partition p for weight dimension w, and the part sizes for
 *     weight dimension w sum to 1.0.
 *
 *  \param  partNums  is the list of part numbers for the weights
 *     which are provided in partWeights.
 *
 *  \param  partWeights  one array for each weight dimension.  
 *     partWeights[w][i] is the weight for part partNums[i]
 *     in weight dimension w.  Other processes may also have
 *     contributions to the weight of partNums[i].
 *
 *  \param result is an array allocated by the caller.
 *     It has space for one imbalance for each weight dimension.
 *     On return, result[w] will have the global imbalance
 *     for weight dimension w.
 *
 *  The weight dimension (partSizes.size(), partWeights.size(),
 *  and result.size()) should be the same on all processes.
 *
 *  For each weight dimension "w", partSizes[w].size() may either
 *  be zero or numGlobalParts.
 *
 *  For each weight dimension "w", partNums.size() should
 *  equal partWeights[w].size(), since the latter
 *  is the weight corresponding to each of the former.
 *
 *  A technicality: It is possible to compute the imbalance of the
 *  current partitioning to a desired partitioning.  In that case
 *  numGlobalParts and partSizes describes the ideal partitioning.
 *  The current partitioning (partNums and partWeights) may have more 
 *  or fewer parts than numGlobalParts.  In either case we compute
 *  a reasonable value for the imbalance in the current partitioning
 *  with respect to the desired partitioning.
 */

template <typename SCALAR>
  void imbalances(const RCP<const Environment> &env, 
    const RCP<const Comm<int> > &comm, 
    size_t numGlobalParts, Array<ArrayView<float> > &partSizes,
    ArrayView<size_t> partNums, 
    Array<ArrayView<SCALAR> > &partWeights, 
    ArrayView<float> result)
{
  // Minimum and maximum actual part numbers

  size_t lmax=0, lmin=numGlobalParts;
  size_t nparts = partNums.size();

  for (size_t p=0; p < nparts; p++){
    if (partNums[p] > lmax) lmax = partNums[p];
    if (partNums[p] < lmin) lmin = partNums[p];
  }

  // Verify that input makes sense.

  size_t wgtDim = result.size();
  int fail=0;

  if (partWeights.size()!=wgtDim || partSizes.size()!=wgtDim)
    fail = 1;

  for (size_t w=0; !fail && (w < wgtDim); w++){
    if (partNums.size() != partWeights[w].size())
      fail = 1;
    if (partSizes[w].size() != 0 && partSizes[w].size() != numGlobalParts)
      fail = 1;
  }

  Tuple<size_t, 5> lval, gval;

  lval[0] = fail;
  lval[1] = wgtDim;
  lval[2] = -wgtDim;
  lval[3] = -lmin;
  lval[4] = lmax;

  reduceAll<int, size_t>(*comm, Teuchos::REDUCE_MAX, 5,
    lval.getRawPtr(), gval.getRawPtr());

  env->globalInputAssertion(__FILE__, __LINE__, "invalid arguments",
    gval[0]==0 && gval[1]==-gval[2], BASIC_ASSERTION, comm);

  // Get global total weight for all parts

  size_t minPart = -gval[3];
  size_t maxPart = gval[4];
  size_t numParts = maxPart - minPart + 1;  // not necessarily numGlobalParts

  Array<SCALAR> myWeights(wgtDim * numParts, 0.0);
  Array<SCALAR> sumWeights(wgtDim * numParts, 0.0);
  
  for (size_t w=0,offset=0; w < wgtDim; w++,offset+=numParts){
    ArrayView<SCALAR> wgt = partWeights[w];
    ArrayView<SCALAR> outWeights = myWeights.view(offset, numParts);

    for (size_t p = 0; p < nparts; p++){
      size_t idx = partNums[p] - minPart;
      outWeights[idx] = wgt[p];
    }
  }

  reduceAll<int, SCALAR>(*comm, Teuchos::REDUCE_SUM, wgtDim * numParts, 
    myWeights.getRawPtr(), sumWeights.getRawPtr());

  SCALAR sizeFactor = 1.0 / numGlobalParts;

  size_t fromPart = 0;
  size_t toPart = numGlobalParts - 1;

  if (minPart > fromPart) fromPart = minPart;
  if (maxPart < toPart)   toPart = maxPart;

  for (size_t w=0, offset=0; w < wgtDim; w++, offset+=numParts){

    ArrayView<SCALAR> inWeights = sumWeights.view(offset, numParts);

    SCALAR total = 0.0, max=0.0;

    for (int p=0; p < numParts; p++)
      total += inWeights[p];

    if (partSizes[w].size() > 0){
      for (size_t p = fromPart; p <= toPart; p++){
        size_t idx = p - minPart;
        SCALAR denom = total*partSizes[w][p]; 
        if (denom > 0){
          SCALAR tmp = inWeights[idx] / denom;
          if (tmp > max) max = tmp;
        }
      }
    }
    else{
      SCALAR denom = total*sizeFactor;
      for (size_t p = fromPart; p <= toPart; p++){
        size_t idx = p - minPart;
        SCALAR tmp = inWeights[idx] / denom;
        if (tmp > max) max = tmp;
      }
    }

    if (max == 0)
      result[w] = -1;
    else
      result[w] = max;
  }
}

/*! \brief Compute the global imbalance.  Number of weights per object is one.
 *
 *  \param  env    library environment information
 *
 *  \param  comm   the communicator
 *
 *  \param  numGlobalParts   the desired global number of parts
 *
 *  \param  partSizes If partSizes.size() is zero, then we assume 
 *     uniform part sizes are desired. Otherwise partSizes has
 *     size numGlobalParts, partSizes[p] is the part size desired for 
 *     partition p, and the part sizes for sum to 1.0.
 *
 *  \param  partNums  is the list of part numbers for the weights
 *     which are provided in partWeights.
 *
 *  \param  partWeights  partWeights[i] is the weight for part partNums[i].
 *     Other processes may also have contributions to the weight of partNums[i].
 *
 *  \param result On return, result will have the global imbalance.
 *
 *  partSizes.size() may either be zero or numGlobalParts.
 *
 *  partNums.size() should equal partWeights.size(), since the latter
 *  is the weight corresponding to each of the former.
 *
 *  A technicality: It is possible to compute the imbalance of the
 *  current partitioning to a desired partitioning.  In that case
 *  numGlobalParts and partSizes describes the ideal partitioning.
 *  The current partitioning (partNums and partWeights) may have more 
 *  or fewer parts than numGlobalParts.  In either case we compute
 *  a reasonable value for the imbalance in the current partitioning
 *  with respect to the desired partitioning.
 */

template <typename SCALAR>
  void imbalances(const RCP<const Environment> &env, 
    const RCP<const Comm<int> > &comm, 
    size_t numGlobalParts, ArrayView<float> partSizes,
    ArrayView<size_t> partNums, ArrayView<SCALAR> partWeights, 
    float &result)
{
  Array<ArrayView<float> > partSizeLists(1, partSizes);
  Array<ArrayView<SCALAR> > partWeightLists(1, partWeights);
  Array<float> resultList(1, 0);

  imbalances(env, comm, numGlobalParts, partSizeLists, partNums,
    partWeightLists, resultList.view(0, 1));

  result = resultList[0];
}

} //namespace Zoltan2
#endif
