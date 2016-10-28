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

#include <Zoltan2_EvaluateBaseClass.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

namespace Zoltan2{

/*! \brief A class that computes and returns quality metrics.
 *  \todo For some problems it will be necessary to build the
 *          Model again in order to compute metrics.  For now
 *          we don't have any problems like that.
    \todo write a unit test for this class
 */

template <typename Adapter>
class EvaluatePartition : public EvaluateBaseClass<Adapter> {

  // defines for metric checks - these are key words used in xml
  #define WEIGHT_PARAMETER_NAME "weight"
  #define NORMED_PARAMETER_NAME "normed"

private:
  using EvaluateBaseClass<Adapter>::getAllMetricsOfType;
  using EvaluateBaseClass<Adapter>::getAllMetrics;

  typedef typename Adapter::base_adapter_t base_adapter_t;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::part_t part_t;
  typedef typename Adapter::scalar_t scalar_t;
  typedef BaseClassMetrics<scalar_t> base_metric_t;

  typedef StridedData<lno_t, scalar_t> input_t;
  part_t numGlobalParts_;           // desired
  part_t targetGlobalParts_;        // actual
  part_t numNonEmpty_;              // of actual

  void sharedConstructor(const Adapter *ia,
                         ParameterList *p,
                         const RCP<const Comm<int> > &problemComm,
                         const PartitioningSolution<Adapter> *soln,
                         const RCP<const GraphModel
                           <typename Adapter::base_adapter_t> > &graphModel);



  /*! \brief Constructor where communicator is Teuchos default, and takes
   *  another parameter whether to evaluate metrics within the constructor or not.
      \param ia the problem input adapter
      \param p the parameter list
      \param problemComm  the problem communicator
      \param soln  the solution
      \param force_evaluate whether to evaluate within to constructor or not.
      \param graphModel the graph model
      The constructor does global communication to compute the metrics .
      The rest of the  methods are local.
   */
  EvaluatePartition(
    const Adapter *ia,
    ParameterList *p,
    const RCP<const Comm<int> > &problemComm,
    const PartitioningSolution<Adapter> *soln,
    bool force_evaluate,
    const RCP<const GraphModel<typename Adapter::base_adapter_t> > &graphModel=
          Teuchos::null):
            numGlobalParts_(0), targetGlobalParts_(0), numNonEmpty_(0), metricsBase_() {
    if (force_evaluate){
      sharedConstructor(ia, p, problemComm, soln, graphModel);
    }

  }
  virtual void calculate_graph_metrics(
      const RCP<const Environment> &_env,
      const RCP<const Comm<int> > &_problemComm,
      const RCP<const GraphModel<typename Adapter::base_adapter_t> > &_graph,
      const ArrayView<const typename Adapter::part_t> &_partArray,
      typename Adapter::part_t &_numGlobalParts,
      ArrayRCP<RCP<BaseClassMetrics<typename Adapter::scalar_t> > > &_metricsBase,
      ArrayRCP<typename Adapter::scalar_t> &_globalSums){
    globalWeightedCutsMessagesByPart <Adapter>(_env,
            _problemComm, _graph, _partArray,
            _numGlobalParts, _metricsBase,
            _globalSums);
  }
public:
  virtual ~EvaluatePartition(){}

  /*! \brief Constructor where communicator is Teuchos default.
      \param ia the problem input adapter
      \param p the parameter list
      \param soln  the solution
      \param graphModel the graph model
      The constructor does global communication to compute the metrics.
      The rest of the  methods are local.
   */
  EvaluatePartition(
    const Adapter *ia, 
    ParameterList *p,
    const PartitioningSolution<Adapter> *soln,
    const RCP<const GraphModel<typename Adapter::base_adapter_t> > &graphModel= 
          Teuchos::null):
    numGlobalParts_(0), targetGlobalParts_(0), numNonEmpty_(0)
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
  EvaluatePartition(
    const Adapter *ia,
    ParameterList *p,
    const RCP<const Comm<int> > &problemComm,
    const PartitioningSolution<Adapter> *soln,
    const RCP<const GraphModel<typename Adapter::base_adapter_t> > &graphModel=
          Teuchos::null):
    numGlobalParts_(0), targetGlobalParts_(0), numNonEmpty_(0)
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
  EvaluatePartition(
    const Adapter *ia, 
    ParameterList *p,
    MPI_Comm comm,
    const PartitioningSolution<Adapter> *soln,
    const RCP<const GraphModel<typename Adapter::base_adapter_t> > &graphModel=
          Teuchos::null):
    numGlobalParts_(0), targetGlobalParts_(0), numNonEmpty_(0)
  {
    RCP<Teuchos::OpaqueWrapper<MPI_Comm> > wrapper =
        Teuchos::opaqueWrapper(comm);
    RCP<const Comm<int> > problemComm =
        rcp<const Comm<int> >(new Teuchos::MpiComm<int>(wrapper));
    sharedConstructor(ia, p, problemComm, soln, graphModel);
  }
#endif

  /*! \brief Return the object count imbalance.
   */
  scalar_t getObjectCountImbalance() const {
    auto metrics = getAllMetricsOfType(IMBALANCE_METRICS_TYPE_NAME);
    if( metrics.size() <= 0 ) {
      throw std::logic_error("getObjectCountImbalance() was called " 
                             "but no metrics data was generated for " + 
                             std::string(IMBALANCE_METRICS_TYPE_NAME) + "." );
    }
    return metrics[0]->getMetricValue("maximum imbalance");
  }

  /*! \brief Return the object normed weight imbalance.
   *  Normed imbalance is only valid if there is at least 2 elements - 
   *  the second one is the normed imbalance.
   *  If we have weights (which start at the second element) the spec is to 
   *  have this return that element.
   */
  scalar_t getNormedImbalance() const{
    auto metrics = getAllMetricsOfType(IMBALANCE_METRICS_TYPE_NAME);
    if( metrics.size() <= 0 ) {
      throw std::logic_error("getNormedImbalance() was called " 
                             "but no metrics data was generated for " + 
                             std::string(IMBALANCE_METRICS_TYPE_NAME) + "." );
    }
    if( metrics.size() <= 1 ) {
      throw std::logic_error("getNormedImbalance() was called " 
                             "but the normed data does not exist." );
    }
    return metrics[1]->getMetricValue("maximum imbalance");
  }

  /*! \brief Return the imbalance for the requested weight.
   */
  scalar_t getWeightImbalance(int weightIndex) const {
    // In this case we could have
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

    // if we have multiple weights (meaning array size if 2 or greater) than 
    // the weights begin on index 2
    // if we have one weight 0 (option 2) then the weights begin on index 1
    auto metrics = getAllMetricsOfType(IMBALANCE_METRICS_TYPE_NAME);
    int weight0IndexStartsAtThisArrayIndex = ( metrics.size() > 2 ) ? 2 : 1;
    int numberOfWeights = metrics.size() - weight0IndexStartsAtThisArrayIndex;
    int indexInArray = weight0IndexStartsAtThisArrayIndex + weightIndex;
    if( metrics.size() <= indexInArray ) {
      throw std::logic_error("getWeightImbalance was called with weight index "+
                             std::to_string(weightIndex) + 
                             " but the maximum weight available for " + 
                             std::string(IMBALANCE_METRICS_TYPE_NAME) + 
                             " is weight " + std::to_string(numberOfWeights-1) +
                             "." );
    }
    return metrics[indexInArray]->getMetricValue("maximum imbalance");
  }

  /*! \brief Return the max cut for the requested weight.
   */
  scalar_t getMaxEdgeCut() const{
    auto graphMetrics = getAllMetricsOfType(GRAPH_METRICS_TYPE_NAME);
    if( graphMetrics.size() < 1 ) {
      throw std::logic_error("getMaxEdgeCut() was called " 
                             "but no metrics data was generated for " + 
                             std::string(GRAPH_METRICS_TYPE_NAME) + "." );
    }
    return graphMetrics[0]->getMetricValue("global maximum");
  }

  /*! \brief getMaxWeightEdgeCuts weighted for the specified index
   */
  scalar_t getMaxWeightEdgeCut(int weightIndex) const{
    auto graphMetrics = getAllMetricsOfType(GRAPH_METRICS_TYPE_NAME);
    int indexInArray = weightIndex + 1; // changed this - it used to start at 0
    if( graphMetrics.size() <= 1 ) {
      throw std::logic_error("getMaxWeightEdgeCut was called with " 
                             "weight index " + std::to_string(weightIndex) + 
                             " but no weights were available for " + 
                             std::string(GRAPH_METRICS_TYPE_NAME) + "." );
    }
    else if( graphMetrics.size() <= indexInArray ) {
      // the size() - 2 is because weight 0 starts at array element 1 
      // (so if the array size is 2, the maximum specified weight index is 
      // weight 0 ( 2-2 = 0 )
      throw std::logic_error("getMaxWeightEdgeCut was called with " 
                             "weight index " + std::to_string(weightIndex) + 
                             " but the maximum weight available for " + 
                             std::string(GRAPH_METRICS_TYPE_NAME) + 
                             " is weight " + 
                             std::to_string(graphMetrics.size() - 2) + "." );
    }
    return graphMetrics[indexInArray]->getMetricValue("global maximum");
  }

  /*! \brief getTotalEdgeCut
   */
  scalar_t getTotalEdgeCut() const{
    auto graphMetrics = getAllMetricsOfType(GRAPH_METRICS_TYPE_NAME);
    if( graphMetrics.size() < 1 ) {
      throw std::logic_error("getTotalEdgeCut() was called but no metrics " 
                             "data was generated for " + 
                             std::string(GRAPH_METRICS_TYPE_NAME) + "." );
    }
    return graphMetrics[0]->getMetricValue("global sum");
  }

  /*! \brief getTotalWeightEdgeCut weighted for the specified index
   */
  scalar_t getTotalWeightEdgeCut(int weightIndex) const{
    auto graphMetrics = getAllMetricsOfType(GRAPH_METRICS_TYPE_NAME);
    int indexInArray = weightIndex + 1; // changed this; it used to start at 0
    if( graphMetrics.size() <= 1 ) { 
      // the size() - 2 is because weight 0 starts at array element 1 (so if 
      // the array size is 2, the maximum specified weight index is 
      // weight 0 ( 2-2 = 0 )
      throw std::logic_error("getTotalWeightEdgeCut was called with " 
                             "weight index " + std::to_string(weightIndex) + 
                             " but no weights were available for " + 
                             std::string(GRAPH_METRICS_TYPE_NAME) + "." );
    }
    else if( graphMetrics.size() <= indexInArray ) {
      throw std::logic_error("getTotalWeightEdgeCut was called with " 
                             "weight index " + std::to_string(weightIndex) + 
                             " but the maximum weight available for " + 
                             std::string(GRAPH_METRICS_TYPE_NAME) + 
                             " is weight " + 
                             std::to_string(graphMetrics.size() - 2) + "." );
    }
    return graphMetrics[indexInArray]->getMetricValue("global sum");
  }

  /*! \brief Print all metrics of type metricType based on the metric object type
   *  Note that parent class currently suppresses this if the list is empty.
   */
  virtual void callStaticPrintMetrics(std::ostream &os,
    ArrayView<RCP<base_metric_t>> metrics, std::string metricType) const {
      if( metricType == GRAPH_METRICS_TYPE_NAME ) {
        Zoltan2::printGraphMetrics<scalar_t, part_t>(os, targetGlobalParts_,
                            numGlobalParts_, metrics);
      }
      else if( metricType == IMBALANCE_METRICS_TYPE_NAME ) {
        Zoltan2::printImbalanceMetrics<scalar_t, part_t>(os, targetGlobalParts_,
                            numGlobalParts_, numNonEmpty_, metrics);
    }
  }

  /*! \brief Return true for any names we accept.
   */
  virtual bool isMetricCheckNameValid(std::string metricCheckName) const {
    return (metricCheckName == WEIGHT_PARAMETER_NAME ||
            metricCheckName == NORMED_PARAMETER_NAME);
  }

  /*! \brief Reads a metric value for bounds checking.
   * Handle any special optional parameters.
   */
  virtual MetricAnalyzerInfo<typename Adapter::scalar_t> getMetricResult(
    const ParameterList & metricCheckParameters, std::string keyWord) const {

    MetricAnalyzerInfo<typename Adapter::scalar_t> result;

    // didn't want to duplicate this value - a weight index should be 0 or
    // larger but it's optional to specify it
    #define UNDEFINED_PARAMETER_INT_INDEX -1

    // Read the weight index parameter and throw if not a good format
    // This is an optional parameter for EvaluatePartition
    int weightIndex = UNDEFINED_PARAMETER_INT_INDEX;
    if( metricCheckParameters.isParameter(WEIGHT_PARAMETER_NAME)) {
      weightIndex = metricCheckParameters.get<int>(WEIGHT_PARAMETER_NAME);
      if( weightIndex < 0 ) {
        throw std::logic_error( "Optional weight index was specified as: " +
          std::to_string(weightIndex) +
          "   Weight index must be 0 or positive." );
      }
    }

    // Read the norm index and throw if not a good format
    // This is an optional parameter for EvaluatePartition
    int normedSetting = UNDEFINED_PARAMETER_INT_INDEX;
    if( metricCheckParameters.isParameter(NORMED_PARAMETER_NAME)) {
      bool bNormSetting = metricCheckParameters.get<bool>(NORMED_PARAMETER_NAME);
      normedSetting = bNormSetting ? 1 : 0;
      if( normedSetting != 0 && normedSetting != 1 ) {
        throw std::logic_error( "Optional normed parameter was specified as: "
          + std::to_string(normedSetting) +
          "   Normed parameter must be true or false." );
      }
    }

    if( weightIndex != UNDEFINED_PARAMETER_INT_INDEX &&
      normedSetting != UNDEFINED_PARAMETER_INT_INDEX ) {
      throw std::logic_error( "Both parameters 'normed' and 'weight' were "
        " specified. They should never appear together." );
    }

    // these define key names which convert to an API call
    #define API_STRING_getWeightImbalance "imbalance"
    #define API_STRING_getTotalEdgeCuts "total edge cuts"
    #define API_STRING_getMaxEdgeCuts "max edge cuts"

    // throw if normed set and weight is set
    if( keyWord != API_STRING_getWeightImbalance &&
      normedSetting != UNDEFINED_PARAMETER_INT_INDEX ) {
      throw std::logic_error( "'normed' was specified but this only has meaning"
       " for the 'imbalance' parameter." );
    }

    // Enforcing parallel usage to the API calls exist in EvaluatePartition
    if (keyWord == API_STRING_getWeightImbalance) {
      if( weightIndex == UNDEFINED_PARAMETER_INT_INDEX ) {
        if( normedSetting == 1 ) {
          result.theValue = getNormedImbalance();
        }
        else {
          result.theValue = getObjectCountImbalance(); // this will be index
        }
      }
      else {
        // this will get the proper index specified
        result.theValue = getWeightImbalance(weightIndex);
      }
    }
    else if (keyWord == API_STRING_getTotalEdgeCuts) {
      if( weightIndex == UNDEFINED_PARAMETER_INT_INDEX ) {
        result.theValue = getTotalEdgeCut();
      }
      else {
        result.theValue = getTotalWeightEdgeCut(weightIndex);
      }
    }
    else if (keyWord == API_STRING_getMaxEdgeCuts) {
      if( weightIndex == UNDEFINED_PARAMETER_INT_INDEX ) {
        result.theValue = getMaxEdgeCut();
      }
      else {
        result.theValue = getMaxWeightEdgeCut(weightIndex);
      }
    }
    else {
      // we have found an invalid key word - throw an error
      throw std::logic_error( "The parameter '" +
        std::string(KEYWORD_PARAMETER_NAME) + "' was specified as '" +
        keyWord + "' which is not understood." );
    }

    result.parameterDescription = keyWord;
    if( weightIndex != UNDEFINED_PARAMETER_INT_INDEX ) {
      result.parameterDescription = result.parameterDescription +
        " (weight: " + std::to_string(weightIndex) + ")";
    }
    else if( normedSetting != UNDEFINED_PARAMETER_INT_INDEX ) {
      // throw above would catch the case where both of these were set
      result.parameterDescription = result.parameterDescription + " (normed: "
        + ( ( normedSetting == 0 ) ? "false" : "true" ) + ")";
    }

    return result;
  }
};

// sharedConstructor
template <typename Adapter>
void EvaluatePartition<Adapter>::sharedConstructor(
  const Adapter *ia, 
  ParameterList *p,
  const RCP<const Comm<int> > &comm,
  const PartitioningSolution<Adapter> *soln,
  const RCP<const GraphModel<typename Adapter::base_adapter_t> > &graphModel
)
{
  RCP<const Comm<int> > problemComm = (comm == Teuchos::null) ?
    DefaultComm<int>::getComm() : comm;

  RCP<Environment> env;

  try{
    env = rcp(new Environment(*p, problemComm));
  }
  Z2_FORWARD_EXCEPTIONS

  env->debug(DETAILED_STATUS, std::string("Entering EvaluatePartition"));
  env->timerStart(MACRO_TIMERS, "Computing metrics");

  // Parts to which objects are assigned.
  size_t numLocalObjects = ia->getLocalNumIDs();
  ArrayRCP<const part_t> parts;

  if (soln) {
    // User provided a partitioning solution; use it.
    parts = arcp(soln->getPartListView(), 0, numLocalObjects, false);
    env->localInputAssertion(__FILE__, __LINE__, "parts not set",
      ((numLocalObjects == 0) || soln->getPartListView()), BASIC_ASSERTION);
  } else {
    // User did not provide a partitioning solution;
    // Use input adapter's partition.
    const part_t *tmp = NULL;
    ia->getPartsView(tmp);
    if (tmp != NULL) 
      parts = arcp(tmp, 0, numLocalObjects, false);
    else {
      // User has not provided input parts in input adapter
      part_t *procs = new part_t[numLocalObjects];
      for (size_t i=0;i<numLocalObjects;i++) procs[i]=problemComm->getRank();
      parts = arcp(procs, 0, numLocalObjects, true);
    }
  }
  ArrayView<const part_t> partArray = parts(0, numLocalObjects);

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

  try{
    imbalanceMetrics<Adapter>(env, problemComm, mcnorm, ia, soln, partArray, 
      graphModel, numGlobalParts_, numNonEmpty_, getAllMetrics());
  }
  Z2_FORWARD_EXCEPTIONS;

  if (soln)
    targetGlobalParts_ = soln->getTargetGlobalNumberOfParts();
  else
    targetGlobalParts_ = problemComm->getSize();

  env->timerStop(MACRO_TIMERS, "Computing metrics");

  const RCP<const base_adapter_t> bia =
    rcp(dynamic_cast<const base_adapter_t *>(ia), false);
  BaseAdapterType inputType = bia->adapterType();

  if (inputType == GraphAdapterType ||
      inputType == MatrixAdapterType ||
      inputType == MeshAdapterType){
    env->timerStart(MACRO_TIMERS, "Computing graph metrics");
    // When we add parameters for which weights to use, we
    // should check those here.  For now we compute graph metrics
    // using all weights.

    // Create a GraphModel based on input data.
    std::bitset<NUM_MODEL_FLAGS> modelFlags;
    RCP<const GraphModel<base_adapter_t> > graph = graphModel;
    if (graph == Teuchos::null) {
      graph = rcp(new GraphModel<base_adapter_t>(bia, env, problemComm,
        modelFlags));
    }

    try {
      ArrayRCP<scalar_t> globalSums; // compute weighted cuts
      globalWeightedCutsByPart<Adapter>(env, problemComm, graph, partArray,
        numGlobalParts_, getAllMetrics(), globalSums);
        
      this->calculate_graph_metrics(env, problemComm, graph, partArray,
        numGlobalParts_, getAllMetrics(), globalSums);
    }
    Z2_FORWARD_EXCEPTIONS;

    env->timerStop(MACRO_TIMERS, "Computing graph metrics");
  }

  env->debug(DETAILED_STATUS, std::string("Exiting EvaluatePartition"));
}

}   // namespace Zoltan2

#endif
