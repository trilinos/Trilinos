// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/*! \file Zoltan2_GraphMetricValuesUtility.hpp
 */

#ifndef ZOLTAN2_GRAPHICMETRICVALUESUTILITY_HPP
#define ZOLTAN2_GRAPHICMETRICVALUESUTILITY_HPP

#include <Zoltan2_ImbalanceMetrics.hpp>
#include <Zoltan2_MetricUtility.hpp>

namespace Zoltan2{

/*! \brief Given the local partitioning, compute the global weighted cuts in each part.
 *
 *   \param env   Environment for error handling
 *   \param comm   communicator
 *   \param graph Graph model
 *   \param part   \c part[i] is the part ID for local object \c i
 *   \param numParts  on return this is the global number of parts.
 *   \param metrics on return points to a list of named GraphMetricValues cuts
 *     that each contains the global max and sum over parts of
 *     the item being measured. The list may contain "cut count", or
 *     "weight 0", "weight 1" and so on in that order.
 *     If uniform weights were given, then only "cut count" appears.
 *     If one set of non-uniform weights were given, then
 *     "weight 0" appear.  Finally, if multiple
 *     weights were given, we have
 *     the individual weights "weight 0", "weight 1", and so on.
 *   \param globalSums If weights are uniform, the globalSums is the
 *      \c numParts totals of global number of cuts in each part.
 *     Suppose the number of weights is \c W.  If
 *     W is 1, then on return this is an array of length \c numParts .
 *     The \c numParts entries are the total weight in each part.
 *     If \c W is greater than one, then the length of this array is
 *     \c W*numParts .
 *     The entries are the sum of the individual weights in each part,
 *     by weight index by part number.  The array is allocated here.
 *
 * globalWeightedCutsByPart() must be called by all processes in \c comm.
 */

template <typename Adapter>
  void globalWeightedCutsByPart(
    const RCP<const Environment> &env,
    const RCP<const Comm<int> > &comm,
    const RCP<const GraphModel<typename Adapter::base_adapter_t> > &graph,
    const ArrayView<const typename Adapter::part_t> &part,
    typename Adapter::part_t &numParts,
    ArrayRCP<RCP<BaseClassMetrics<typename Adapter::scalar_t> > > &metrics,
    ArrayRCP<typename Adapter::scalar_t> &globalSums)
{
  env->debug(DETAILED_STATUS, "Entering globalWeightedCutsByPart");
  //////////////////////////////////////////////////////////
  // Initialize return values

  numParts = 0;

  int ewgtDim = graph->getNumWeightsPerEdge();

  int numMetrics = 1;                   // "cut count"
  if (ewgtDim) numMetrics += ewgtDim;   // "weight n"

  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::node_t node_t;
  typedef typename Adapter::part_t part_t;
  typedef StridedData<lno_t, scalar_t> input_t;

  typedef GraphMetrics<scalar_t> mv_t;
  typedef Tpetra::CrsMatrix<part_t,lno_t,gno_t,node_t>  sparse_matrix_type;
  typedef Tpetra::Vector<part_t,lno_t,gno_t,node_t>     vector_t;
  typedef Tpetra::Map<lno_t, gno_t, node_t>                map_type;
  typedef Tpetra::global_size_t GST;
  const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

  using Teuchos::as;

  // add some more metrics to the array
  typedef typename ArrayRCP<RCP<BaseClassMetrics<typename Adapter::scalar_t> > >::size_type array_size_type;
  metrics.resize( metrics.size() + numMetrics );
  for( array_size_type n = metrics.size() - numMetrics; n < metrics.size(); ++n )
  {
	  mv_t * newMetric = new mv_t;									// allocate the new memory
	  env->localMemoryAssertion(__FILE__,__LINE__,1,newMetric);		// check errors
	  metrics[n] = rcp( newMetric); 				// create the new members
  }
  array_size_type next = metrics.size() - numMetrics; // MDM - this is most likely temporary to preserve the format here - we are now filling a larger array so we may not have started at 0


  //////////////////////////////////////////////////////////
  // Figure out the global number of parts in use.
  // Verify number of vertex weights is the same everywhere.

  lno_t localNumObj = part.size();
  part_t localNum[2], globalNum[2];
  localNum[0] = static_cast<part_t>(ewgtDim);
  localNum[1] = 0;

  for (lno_t i=0; i < localNumObj; i++)
    if (part[i] > localNum[1]) localNum[1] = part[i];

  try{
    reduceAll<int, part_t>(*comm, Teuchos::REDUCE_MAX, 2,
      localNum, globalNum);
  }
  Z2_THROW_OUTSIDE_ERROR(*env)

  env->globalBugAssertion(__FILE__,__LINE__,
    "inconsistent number of edge weights",
    globalNum[0] == localNum[0], DEBUG_MODE_ASSERTION, comm);

  part_t nparts = globalNum[1] + 1;

  part_t globalSumSize = nparts * numMetrics;
  scalar_t * sumBuf = new scalar_t [globalSumSize];
  env->localMemoryAssertion(__FILE__, __LINE__, globalSumSize, sumBuf);
  globalSums = arcp(sumBuf, 0, globalSumSize);

  //////////////////////////////////////////////////////////
  // Calculate the local totals by part.

  scalar_t *localBuf = new scalar_t [globalSumSize];
  env->localMemoryAssertion(__FILE__,__LINE__,globalSumSize,localBuf);
  memset(localBuf, 0, sizeof(scalar_t) * globalSumSize);

  scalar_t *cut = localBuf;              // # of cuts

  ArrayView<const gno_t> Ids;
  ArrayView<input_t> vwgts;
  //size_t nv =
  graph->getVertexList(Ids, vwgts);

  ArrayView<const gno_t> edgeIds;
  ArrayView<const lno_t> offsets;
  ArrayView<input_t> wgts;
  //size_t numLocalEdges =
  graph->getEdgeList(edgeIds, offsets, wgts);
  // **************************************************************************
  // *************************** BUILD MAP FOR ADJS ***************************
  // **************************************************************************

  RCP<const map_type> vertexMapG;

  // Build a list of the global vertex ids...
  gno_t min = std::numeric_limits<gno_t>::max();
  size_t maxcols = 0;
  for (lno_t i = 0; i < localNumObj; ++i) {
    if (Ids[i] < min) min = Ids[i];
    size_t ncols = offsets[i+1] - offsets[i];
    if (ncols > maxcols) maxcols = ncols;
  }

  gno_t gmin;
  Teuchos::reduceAll<int, gno_t>(*comm,Teuchos::REDUCE_MIN,1,&min,&gmin);

  //Generate Map for vertex
  vertexMapG = rcp(new map_type(INVALID, Ids, gmin, comm));

  // **************************************************************************
  // ************************** BUILD GRAPH FOR ADJS **************************
  // **************************************************************************

  RCP<sparse_matrix_type> adjsMatrix;

  // Construct Tpetra::CrsGraph objects.
  adjsMatrix = rcp (new sparse_matrix_type (vertexMapG, 0));

  Array<part_t> justOneA(maxcols, 1);

  for (lno_t localElement=0; localElement<localNumObj; ++localElement){
    // Insert all columns for global row Ids[localElement]
    size_t ncols = offsets[localElement+1] - offsets[localElement];
    adjsMatrix->insertGlobalValues(Ids[localElement],
                                   edgeIds(offsets[localElement], ncols),
                                   justOneA(0, ncols));
  }

  //Fill-complete adjs Graph
  adjsMatrix->fillComplete ();

  // Compute part
  RCP<vector_t> scaleVec = Teuchos::rcp( new vector_t(vertexMapG,false) );
  for (lno_t localElement=0; localElement<localNumObj; ++localElement) {
    scaleVec->replaceLocalValue(localElement,part[localElement]);
  }

  // Postmultiply adjsMatrix by part
  adjsMatrix->rightScale(*scaleVec);
  Array<gno_t> Indices;
  Array<part_t> Values;

  for (lno_t i=0; i < localNumObj; i++) {
    const gno_t globalRow = Ids[i];
    size_t NumEntries = adjsMatrix->getNumEntriesInGlobalRow (globalRow);
    Indices.resize (NumEntries);
    Values.resize (NumEntries);
    adjsMatrix->getGlobalRowCopy (globalRow,Indices(),Values(),NumEntries);

    for (size_t j=0; j < NumEntries; j++)
      if (part[i] != Values[j])
	cut[part[i]]++;
  }

  if (numMetrics > 1) {

    scalar_t *wgt = localBuf + nparts; // weight 0

    // This code assumes the solution has the part ordered the
    // same way as the user input.  (Bug 5891 is resolved.)
    for (int edim = 0; edim < ewgtDim; edim++){
      for (lno_t i=0; i < localNumObj; i++) {
	const gno_t globalRow = Ids[i];
	size_t NumEntries = adjsMatrix->getNumEntriesInGlobalRow (globalRow);
	Indices.resize (NumEntries);
	Values.resize (NumEntries);
	adjsMatrix->getGlobalRowCopy (globalRow,Indices(),Values(),NumEntries);

	for (size_t j=0; j < NumEntries; j++)
	  if (part[i] != Values[j])
	    wgt[part[i]] += wgts[edim][offsets[i] + j];
      }
      wgt += nparts;         // individual weights
    }
  }

  //////////////////////////////////////////////////////////
  // Obtain global totals by part.

  try{
    reduceAll<int, scalar_t>(*comm, Teuchos::REDUCE_SUM, globalSumSize,
      localBuf, sumBuf);
  }
  Z2_THROW_OUTSIDE_ERROR(*env);

  delete [] localBuf;

  //////////////////////////////////////////////////////////
  // Global max and sum over all parts

  cut = sumBuf;                     // # of cuts
  scalar_t max=0, sum=0;

  ArrayView<scalar_t> cutVec(cut, nparts);
  getStridedStats<scalar_t>(cutVec, 1, 0, max, sum);

  metrics[next]->setName("cut count");
  metrics[next]->setMetricValue("global maximum", max);
  metrics[next]->setMetricValue("global sum", sum);

  next++;

  if (numMetrics > 1){
    scalar_t *wgt = sumBuf + nparts;        // weight 0

    for (int edim=0; edim < ewgtDim; edim++){
      ArrayView<scalar_t> fromVec(wgt, nparts);
      getStridedStats<scalar_t>(fromVec, 1, 0, max, sum);

      std::ostringstream oss;
      oss << "weight " << edim;

      metrics[next]->setName(oss.str());
      metrics[next]->setMetricValue("global maximum", max);
      metrics[next]->setMetricValue("global sum", sum);

      next++;
      wgt += nparts;       // individual weights
    }
  }

  numParts = nparts;

  env->debug(DETAILED_STATUS, "Exiting globalWeightedCutsByPart");
}

/*! \brief Print out a header and the values for a list of graph metrics.
 */
template <typename scalar_t, typename part_t>
  void printMetrics( std::ostream &os,
    part_t targetNumParts, part_t numParts,
    const ArrayRCP<RCP<BaseClassMetrics<scalar_t>>> &infoList)
{
  os << "NUMBER OF PARTS IS " << numParts;
  os << std::endl;
  if (targetNumParts != numParts)
    os << "TARGET NUMBER OF PARTS IS " << targetNumParts << std::endl;

  std::string unset(METRICS_UNSET_STRING);

  GraphMetrics<scalar_t>::printHeader(os);

  for (int i=0; i < infoList.size(); i++)
    if (infoList[i]->getName() != unset)
      infoList[i]->printLine(os);

  os << std::endl;
}

/*! \brief Print out a header and the values for a single metric.
 */
template <typename scalar_t, typename part_t>
  void printMetrics( std::ostream &os,
    part_t targetNumParts, part_t numParts,
    RCP<BaseClassMetrics<scalar_t>> metricValue)
{
  ArrayRCP<RCP<BaseClassMetrics<scalar_t> > > infoList = arcp<RCP<BaseClassMetrics<scalar_t>>>( 1 );	// new
  infoList[0] = metricValue;
  printMetrics( os, targetNumParts, numParts, infoList);
}

} //namespace Zoltan2

#endif
