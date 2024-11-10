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

#ifndef ZOLTAN2_EVALUATEMAPPING_HPP
#define ZOLTAN2_EVALUATEMAPPING_HPP

#include <Zoltan2_MappingSolution.hpp>
#include <Zoltan2_EvaluatePartition.hpp>
#include <Zoltan2_MachineRepresentation.hpp>

namespace Zoltan2{

/*! \brief A class that computes and returns quality metrics.
 *  \todo For some problems it will be necessary to build the
 *          Model again in order to compute metrics.  For now
 *          we don't have any problems like that.
    \todo write a unit test for this class
 */

template<typename Adapter,
         typename MachineRep =   // Default MachineRep type
                  MachineRepresentation<typename Adapter::scalar_t,
                                        typename Adapter::part_t> >
class EvaluateMapping : public EvaluatePartition<Adapter> {


  const RCP <const MachineRep> machine;
public:

  /*! \brief Constructor where communicator is Teuchos default.
      \param ia the problem input adapter
      \param p the parameter list
      \param soln  the solution
      \param graphModel the graph model
      The constructor does global communication to compute the metrics.
      The rest of the  methods are local.
   */
  EvaluateMapping(
    const Adapter *ia, 
    ParameterList *p,
    const RCP<const Comm<int> > &ucomm_,
    const MappingSolution<Adapter> *soln,
    const MachineRep *machine_ ,
    const RCP<const GraphModel<typename Adapter::base_adapter_t> > &graphModel= Teuchos::null):
      EvaluatePartition<Adapter>(ia, p, ucomm_, soln, false, graphModel),
      machine(Teuchos::rcp(machine_, false)){
    this->sharedConstructor(ia, p, ucomm_, soln, graphModel);
  }

  virtual ~EvaluateMapping(){}

protected:
  virtual void calculate_graph_metrics(
      const RCP<const Environment> &_env,
      const RCP<const Comm<int> > &_problemComm,
      const RCP<const GraphModel<typename Adapter::base_adapter_t> > &_graph,
      const ArrayView<const typename Adapter::part_t> &_partArray,
      typename Adapter::part_t &_numGlobalParts,
      ArrayRCP<RCP<BaseClassMetrics<typename Adapter::scalar_t> > > &_metricsBase,
      ArrayRCP<typename Adapter::scalar_t> &_globalSums) {
        globalWeightedByPart <Adapter,MachineRep>(_env,
          _problemComm, _graph, _partArray, _numGlobalParts, _metricsBase,
          _globalSums, true, this->machine);
      }
};

}   // namespace Zoltan2

#endif
