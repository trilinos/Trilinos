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

/*! \file Zoltan2_EvaluatePartition.hpp
 *  \brief Defines the EvaluatePartition class.
 */

#ifndef ZOLTAN2_EVALUATEPARTITION_HPP
#define ZOLTAN2_EVALUATEPARTITION_HPP

#include <Zoltan2_Metric.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

namespace Zoltan2{

/*! \brief A class that computes and returns quality metrics.
 *  \todo For some problems it will be necessary to build the
 *          Model again in order to compute metrics.  For now
 *          we don't have any problems like that.
    \todo write a unit test for this class
 */

  template <typename Adapter>
  class EvaluatePartition {

private:

  typedef typename Adapter::base_adapter_t base_adapter_t;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::part_t part_t;
  typedef typename Adapter::scalar_t scalar_t;

  typedef StridedData<lno_t, scalar_t> input_t;

  part_t numGlobalParts_;           // desired
  part_t targetGlobalParts_;        // actual
  part_t numNonEmpty_;              // of actual

  typedef BaseClassMetrics<scalar_t> base_metric_type;
  typedef ArrayRCP<RCP<base_metric_type>> base_metric_array_type;

  base_metric_array_type metricsBase_;

  void sharedConstructor(const Adapter *ia,
			 ParameterList *p,
			 const RCP<const Comm<int> > &problemComm,
			 const PartitioningSolution<Adapter> *soln,
			 const RCP<const GraphModel
			 <typename Adapter::base_adapter_t> > &graphModel);

public:

  /*! \brief Constructor where communicator is Teuchos default.
      \param ia the problem input adapter
      \param p the parameter list
      \param soln  the solution
      \param graphModel the graph model
      The constructor does global communication to compute the metrics.
      The rest of the  methods are local.
   */
    EvaluatePartition(const Adapter *ia, 
    ParameterList *p,
    const PartitioningSolution<Adapter> *soln,
    const RCP<const GraphModel<typename Adapter::base_adapter_t> > &graphModel=
		    Teuchos::null):
    numGlobalParts_(0), targetGlobalParts_(0), numNonEmpty_(0), metricsBase_()
    {
      RCP<const Comm<int> > problemComm = DefaultComm<int>::getComm();
      sharedConstructor(ia, p, problemComm, soln, graphModel);
    }


  /*! \brief Constructor where Teuchos communicator is specified
      \param ia the problem input adapter
      \param p the parameter list
      \param problemComm  the problem communicator
      \param soln  the solution
      \param graphModel the graph model
      The constructor does global communication to compute the metrics.
      The rest of the  methods are local.
   */
  EvaluatePartition(const Adapter *ia,
    ParameterList *p,
    const RCP<const Comm<int> > &problemComm,
    const PartitioningSolution<Adapter> *soln,
    const RCP<const GraphModel<typename Adapter::base_adapter_t> > &graphModel=
		    Teuchos::null):
    numGlobalParts_(0), targetGlobalParts_(0), numNonEmpty_(0), metricsBase_()
    {
      sharedConstructor(ia, p, problemComm, soln, graphModel);
    }

#ifdef HAVE_ZOLTAN2_MPI
  /*! \brief Constructor for MPI builds
      \param ia the problem input adapter
      \param p the parameter list
      \param comm  the problem communicator
      \param soln  the solution
      \param graphModel the graph model
      The constructor does global communication to compute the metrics.
      The rest of the  methods are local.
   */
    EvaluatePartition(const Adapter *ia, 
    ParameterList *p,
    MPI_Comm comm,
    const PartitioningSolution<Adapter> *soln,
    const RCP<const GraphModel<typename Adapter::base_adapter_t> > &graphModel=
		    Teuchos::null):
    numGlobalParts_(0), targetGlobalParts_(0), numNonEmpty_(0), metricsBase_()
    {
      RCP<Teuchos::OpaqueWrapper<MPI_Comm> > wrapper =
	Teuchos::opaqueWrapper(comm);
      RCP<const Comm<int> > problemComm =
	rcp<const Comm<int> >(new Teuchos::MpiComm<int>(wrapper));
      sharedConstructor(ia, p, problemComm, soln, graphModel);
    }
#endif

  /*! \brief Return the full metric list which will be mixed types.
   */
  base_metric_array_type getAllMetrics() const {
      return metricsBase_;
  }

  /*! \get the unique class names.
   */
  std::set<std::string> getMetricTypes() const {
	  std::set<std::string> metricTypes;
	  for( int n = 0; n < metricsBase_.size(); ++n )
		  metricTypes.insert( metricsBase_[n]->getMetricType() );
	  return metricTypes;
  }

  /*! \brief Return the  metric list for types matching the given metric type.
   */
  base_metric_array_type getAllMetricsOfType(std::string metricType) const
  {
	int counter = 0;
	for( int n = 0; n < metricsBase_.size(); ++n ) {
      if( metricsBase_[n]->getMetricType() == metricType )
      ++counter;
	}

	base_metric_array_type newArray = arcp<RCP<base_metric_type> >( counter );	// new array allocated

	int insertIndex = 0;
	for( int n = 0; n < metricsBase_.size(); ++n ) {
	  if( metricsBase_[n]->getMetricType() == metricType ) {
	    newArray[insertIndex] = metricsBase_[n];
		++insertIndex;
	  }
	}

	return newArray;
  }

  /*! \brief Return the object count imbalance.
   */
  scalar_t getObjectCountImbalance() const {
	base_metric_array_type metrics = getAllMetricsOfType(IMBALANCE_METRICS_TYPE_NAME);

	if( metrics.size() <= 0 ) {
	    throw std::logic_error( "getObjectCountImbalance() was called but no metrics data was generated for " + std::string(IMBALANCE_METRICS_TYPE_NAME) + "." );
	}
	return metrics[0]->getMetricValue("maximum imbalance");
  }

  /*! \brief Return the object normed weight imbalance.
   *  Normed imbalance is only valid if there is at least 2 elements - the second one is the normed imbalance.
   *  If we have weights (which start at the second element) the spec is to have this return that element.
   */
  scalar_t getNormedImbalance() const{
	base_metric_array_type metrics = getAllMetricsOfType(IMBALANCE_METRICS_TYPE_NAME);

	if( metrics.size() <= 0 ) {
	    throw std::logic_error( "getNormedImbalance() was called but no metrics data was generated for " + std::string(IMBALANCE_METRICS_TYPE_NAME) + "." );
	}
	if( metrics.size() <= 1 ) {
		 throw std::logic_error( "getNormedImbalance() was called but the normed data does not exist for " + std::to_string(metrics.size()) + "." );
	}

    return metrics[1]->getMetricValue("maximum imbalance");
  }

  /*! \brief Return the imbalance for the requested weight.
   */
  scalar_t getWeightImbalance(int weightIndex) const {
	// In this case we could have:
		// Option 1
			// object count
		// Option 2
			// object count
			// weight 0
		// Option 3
			// object count
			// normed imbalance
			// weight 0
			// weight 1

	// if we have multiple weights (meaning array size if 2 or greater) than the weights begin on index 2
	// if we have one weight 0 (option 2) then the weights begin on index 1
	base_metric_array_type metrics = getAllMetricsOfType(IMBALANCE_METRICS_TYPE_NAME);
	int weight0IndexStartsAtThisArrayIndex = ( metrics.size() > 2 ) ? 2 : 1;
	int numberOfWeights = metrics.size() - weight0IndexStartsAtThisArrayIndex;
	int useArayIndex = weight0IndexStartsAtThisArrayIndex + weightIndex;
	if( metrics.size() < useArayIndex ) {
   	    throw std::logic_error( "getWeightImbalance(int weightIndex) was called with weight index " + std::to_string(weightIndex) + " but the maximum weight available for " + std::string(IMBALANCE_METRICS_TYPE_NAME) + " is weight " + std::to_string(numberOfWeights-1) + "." );
    }
    return metrics[useArayIndex]->getMetricValue("maximum imbalance");
  }

  /*! \brief Return the max cut for the requested weight.
   */
  scalar_t getMaxEdgeCut() const{
	base_metric_array_type graphMetrics = getAllMetricsOfType(GRAPH_METRICS_TYPE_NAME);
	if( graphMetrics.size() < 1 ) {
	    throw std::logic_error( "getMaxEdgeCut() was called but no metrics data was generated for " + std::string(GRAPH_METRICS_TYPE_NAME) + "." );
	}
	return graphMetrics[0]->getMetricValue("global maximum");
  }

  /*! \brief getMaxWeightEdgeCuts weighted for the specified index
   */
  scalar_t getMaxWeightEdgeCut(int weightIndex) const{
	base_metric_array_type graphMetrics = getAllMetricsOfType(GRAPH_METRICS_TYPE_NAME);
	int indexInArray = weightIndex + 1; // this was changed - it used to start at 0
	if( graphMetrics.size() <= 1 ) {
	    throw std::logic_error( "getMaxWeightEdgeCut(int weightIndex) was called with weight index " + std::to_string(weightIndex) + " but no weights were available for " + std::string(GRAPH_METRICS_TYPE_NAME) + "." );
	}
	else if( graphMetrics.size() < indexInArray ) { // the size() - 2 is because weight 0 starts at array element 1 (so if the array size is 2, the maximum specified weight index is weight 0 ( 2-2 = 0 )
	    throw std::logic_error( "getMaxWeightEdgeCut(int weightIndex) was called with weight index " + std::to_string(weightIndex) + " but the maximum weight available for " + std::string(GRAPH_METRICS_TYPE_NAME) + " is weight " + std::to_string(graphMetrics.size() - 2) + "." );
	}
    return graphMetrics[weightIndex]->getMetricValue("global maximum");
  }

  /*! \brief getTotalEdgeCut
   */
  scalar_t getTotalEdgeCut() const{
	base_metric_array_type graphMetrics = getAllMetricsOfType(GRAPH_METRICS_TYPE_NAME);
	if( graphMetrics.size() < 1 ) {
	    throw std::logic_error( "getTotalEdgeCut() was called but no metrics data was generated for " + std::string(GRAPH_METRICS_TYPE_NAME) + "." );
	}
	return graphMetrics[0]->getMetricValue("global sum");
  }

  /*! \brief getTotalWeightEdgeCut weighted for the specified index
   */
  scalar_t getTotalWeightEdgeCut(int weightIndex) const{
	base_metric_array_type graphMetrics = getAllMetricsOfType(GRAPH_METRICS_TYPE_NAME);
	int indexInArray = weightIndex + 1; // this was changed - it used to start at 0
	if( graphMetrics.size() <= 1 ) { // the size() - 2 is because weight 0 starts at array element 1 (so if the array size is 2, the maximum specified weight index is weight 0 ( 2-2 = 0 )
	    throw std::logic_error( "getTotalWeightEdgeCut(int weightIndex) was called with weight index " + std::to_string(weightIndex) + " but no weights were available for " + std::string(GRAPH_METRICS_TYPE_NAME) + "." );
	}
	else if( graphMetrics.size() < indexInArray ) {
	    throw std::logic_error( "getTotalWeightEdgeCut(int weightIndex) was called with weight index " + std::to_string(weightIndex) + " but the maximum weight available for " + std::string(GRAPH_METRICS_TYPE_NAME) + " is weight " + std::to_string(graphMetrics.size() - 2) + "." );
	}
    return graphMetrics[weightIndex]->getMetricValue("global sum");
  }

  /*! \brief Print all the metrics based on the  metric object type
   */
  void printMetrics(std::ostream &os) const {
	  std::set<std::string> metric_types = getMetricTypes();
	  for( auto metricType : metric_types )
		  printMetrics( os, metricType );
  }

  /*! \brief Print all the metrics of type metricType based on the metric object type
   */
  void printMetrics(std::ostream &os, std::string metricType) const {
	  /*! \Might need changing.
	   *  \The issue here is that they each have unique parameters to pass.
	   */
	  base_metric_array_type metrics = getAllMetricsOfType( metricType );
	  if( metricType == GRAPH_METRICS_TYPE_NAME )
	    Zoltan2::printMetrics<scalar_t, part_t>(os, targetGlobalParts_, numGlobalParts_, metrics);
	  else if( metricType == IMBALANCE_METRICS_TYPE_NAME )
		Zoltan2::printMetrics<scalar_t, part_t>(os, targetGlobalParts_, numGlobalParts_, numNonEmpty_, metrics);
  }
};

  // sharedConstructor
  template <typename Adapter>
  void EvaluatePartition<Adapter>::sharedConstructor(
  const Adapter *ia, 
  ParameterList *p,
  const RCP<const Comm<int> > &comm,
  const PartitioningSolution<Adapter> *soln,
  const RCP<const GraphModel<typename Adapter::base_adapter_t> > &graphModel)
{
  RCP<const Comm<int> > problemComm;
  if (comm == Teuchos::null) {
    problemComm = DefaultComm<int>::getComm();//communicator is Teuchos default
  } else {
    problemComm = comm;
  }
  RCP<Environment> env = rcp(new Environment(*p, problemComm));
  env->debug(DETAILED_STATUS, std::string("Entering EvaluatePartition"));
  env->timerStart(MACRO_TIMERS, "Computing metrics");

  // When we add parameters for which weights to use, we
  // should check those here.  For now we compute metrics
  // using all weights.

  multiCriteriaNorm mcnorm = normBalanceTotalMaximum;
  const Teuchos::ParameterEntry *pe = p->getEntryPtr("partitioning_objective");
  if (pe){
    std::string strChoice = pe->getValue<std::string>(&strChoice);
    if (strChoice == std::string("multicriteria_minimize_total_weight"))
      mcnorm = normMinimizeTotalWeight;
    else if (strChoice == std::string("multicriteria_minimize_maximum_weight"))
      mcnorm = normMinimizeMaximumWeight;
  }

  const RCP<const base_adapter_t> bia =
    rcp(dynamic_cast<const base_adapter_t *>(ia), false);

  try{
    objectMetrics<Adapter>(env, problemComm, mcnorm, ia, soln, graphModel,
			   numGlobalParts_, numNonEmpty_,metricsBase_);
  }


  Z2_FORWARD_EXCEPTIONS;

  if (soln)
  targetGlobalParts_ = soln->getTargetGlobalNumberOfParts();
  else
    targetGlobalParts_ = problemComm->getSize();

  env->timerStop(MACRO_TIMERS, "Computing metrics");

  BaseAdapterType inputType = bia->adapterType();

  if (inputType == GraphAdapterType ||
      inputType == MatrixAdapterType ||
      inputType == MeshAdapterType){
    env->timerStart(MACRO_TIMERS, "Computing graph metrics");
    // When we add parameters for which weights to use, we
    // should check those here.  For now we compute graph metrics
    // using all weights.

    std::bitset<NUM_MODEL_FLAGS> modelFlags;

    // Create a GraphModel based on input data.

    RCP<GraphModel<base_adapter_t> > graph;
    if (graphModel == Teuchos::null)
      graph = rcp(new GraphModel<base_adapter_t>(bia, env, problemComm,
						 modelFlags));

    // Local number of objects.

    size_t numLocalObjects = bia->getLocalNumIDs();

    // Parts to which objects are assigned.

    const part_t *parts;
    if (soln) {
      // User provided a partitioning solution; use it.
      parts = soln->getPartListView();
      env->localInputAssertion(__FILE__, __LINE__, "parts not set",
        ((numLocalObjects == 0) || parts), BASIC_ASSERTION);
    } else {
      // User did not provide a partitioning solution;
      // Use input adapter partition.

      parts = NULL;
      bia->getPartsView(parts);
      if (parts == NULL) {
	// User has not provided input parts in input adapter
	part_t *procs = new part_t [numLocalObjects];
	for (size_t i=0;i<numLocalObjects;i++) procs[i]=problemComm->getRank();
	parts = procs;
      }
    }
    ArrayView<const part_t> partArray(parts, numLocalObjects);

    ArrayRCP<scalar_t> globalSums;

    if (graphModel == Teuchos::null) {
      try{
    	  globalWeightedCutsByPart<Adapter>(env,
					  problemComm, graph, partArray,
					  numGlobalParts_, metricsBase_,
					  globalSums);
      }
      Z2_FORWARD_EXCEPTIONS;
    } else {
      try{
    	  globalWeightedCutsByPart<Adapter>(env,
					  problemComm, graphModel, partArray,
					  numGlobalParts_, metricsBase_,
					  globalSums);
      }
      Z2_FORWARD_EXCEPTIONS;
    }

    env->timerStop(MACRO_TIMERS, "Computing graph metrics");
  }
  env->debug(DETAILED_STATUS,
	     std::string("Exiting EvaluatePartition"));
}

}   // namespace Zoltan2

#endif
