// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Zoltan2_EvaluatePartition.hpp
 *  \brief Defines the EvaluatePartition class.
 */

#ifndef ZOLTAN2_EVALUATEPARTITION_HPP
#define ZOLTAN2_EVALUATEPARTITION_HPP

#include <Zoltan2_GraphMetrics.hpp>
#include <Zoltan2_GraphMetricsUtility.hpp>
#include <Zoltan2_ImbalanceMetrics.hpp>
#include <Zoltan2_ImbalanceMetricsUtility.hpp>
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
  typedef ArrayRCP<RCP<base_metric_type> > base_metric_array_type;
  base_metric_array_type metricsBase_;

protected:
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
      ArrayRCP<typename Adapter::scalar_t> &_globalSums) {
        globalWeightedByPart <Adapter>(_env, _problemComm, _graph,
          _partArray, _numGlobalParts, _metricsBase, _globalSums);
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
    numGlobalParts_(0), targetGlobalParts_(0), numNonEmpty_(0), metricsBase_()
  {
    Teuchos::RCP<const Comm<int> > problemComm = Tpetra::getDefaultComm();
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
  EvaluatePartition(
    const Adapter *ia, 
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

  /*! \brief Return the metric list for types matching the given metric type.
   */

  // TO DO - note that with current status it probably makes more sense to
  // break up metricsBase_ into imbalanceMetrics_ and graphMetrics_.
  // So instead of mixing them just keep two arrays and eliminate this function.
  // That will clean up several places in this class.
  ArrayView<RCP<base_metric_type>> getAllMetricsOfType(
    std::string metricType) const {
    // find the beginning and the end of the contiguous block
    // the list is an ArrayRCP and must preserve any ordering
    int beginIndex = -1;
    int sizeOfArrayView = 0;
    for(auto n = 0; n < metricsBase_.size(); ++n) {
      if( metricsBase_[n]->getMetricType() == metricType ) {
        if (beginIndex == -1) {
          beginIndex = int(n);
        }
        ++sizeOfArrayView;
      }
    }
    if (sizeOfArrayView == 0) {
      return ArrayView<RCP<base_metric_type> >(); // empty array view
    }
    return metricsBase_.view(beginIndex, sizeOfArrayView);
  }

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

  /*! \brief getTotalMessages
   */
  scalar_t getTotalMessages() const{
    auto graphMetrics = getAllMetricsOfType(GRAPH_METRICS_TYPE_NAME);
    if( graphMetrics.size() < 1 ) {
      throw std::logic_error("getTotalMessages() was called but no metrics "
                             "data was generated for " +
                             std::string(GRAPH_METRICS_TYPE_NAME) + "." );
    }
    // TODO: Would be better to avoid hard coding the array access to [1]
    return graphMetrics[1]->getMetricValue("global sum");
  }


  /*! \brief getMaxMessages
   */
  scalar_t getMaxMessages() const{
    auto graphMetrics = getAllMetricsOfType(GRAPH_METRICS_TYPE_NAME);
    if( graphMetrics.size() < 1 ) {
      throw std::logic_error("getMaxMessages() was called but no metrics "
                             "data was generated for " +
                             std::string(GRAPH_METRICS_TYPE_NAME) + "." );
    }
    // TODO: Would be better to avoid hard coding the array access to [1]
    return graphMetrics[1]->getMetricValue("global maximum");
  }

  /*! \brief Print all metrics
   */
  void printMetrics(std::ostream &os) const {
    // this could be a critical decision - do we want a blank table with
    // headers when the list is empty - for debugging that is probably better
    // but it's very messy to have lots of empty tables in the logs
    ArrayView<RCP<base_metric_type>> graphMetrics =
      getAllMetricsOfType(GRAPH_METRICS_TYPE_NAME);
    if (graphMetrics.size() != 0) {
        Zoltan2::printGraphMetrics<scalar_t, part_t>(os, targetGlobalParts_,
          numGlobalParts_, graphMetrics);
    }

    ArrayView<RCP<base_metric_type>> imbalanceMetrics =
      getAllMetricsOfType(IMBALANCE_METRICS_TYPE_NAME);
    if (imbalanceMetrics.size() != 0) {
        Zoltan2::printImbalanceMetrics<scalar_t, part_t>(os, targetGlobalParts_,
          numGlobalParts_, numNonEmpty_, imbalanceMetrics);
    }
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
  RCP<const Comm<int> > problemComm;
  if (comm == Teuchos::null) {
    problemComm = Tpetra::getDefaultComm();
  } else {
    problemComm = comm;
  }

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

  const RCP<const base_adapter_t> bia =
    rcp(dynamic_cast<const base_adapter_t *>(ia), false);

  try{
    imbalanceMetrics<Adapter>(env, problemComm, mcnorm, ia, soln, partArray, 
                              graphModel,
                              numGlobalParts_, numNonEmpty_, metricsBase_);
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

    RCP<const GraphModel<base_adapter_t> > graph = graphModel;
    if (graphModel == Teuchos::null) {
      graph = rcp(new GraphModel<base_adapter_t>(bia, env, problemComm,
        modelFlags));
    }

    // compute weighted cuts
    ArrayRCP<scalar_t> globalSums;
    try {
      this->calculate_graph_metrics(env, problemComm, graph, partArray,
        numGlobalParts_, metricsBase_, globalSums);
    }
    Z2_FORWARD_EXCEPTIONS;

    env->timerStop(MACRO_TIMERS, "Computing graph metrics");
  }

  env->debug(DETAILED_STATUS, std::string("Exiting EvaluatePartition"));
}

}   // namespace Zoltan2

#endif
