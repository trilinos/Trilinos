// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_DAG_MANAGER_DEF_HPP
#define PHX_DAG_MANAGER_DEF_HPP

#include <cstdlib> // for std::getenv
#include <cstddef>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <utility>
#include <typeinfo>
#include "Teuchos_Assert.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator.hpp"
#include "Phalanx_Exceptions.hpp"
#include "Phalanx_FieldTag_STL_Functors.hpp"

#ifdef PHX_ENABLE_KOKKOS_AMT
#include "Kokkos_TaskScheduler.hpp"
#endif

//=======================================================================
template<typename Traits>
PHX::DagManager<Traits>::
DagManager(const std::string& evaluation_type_name) :
  graphviz_filename_for_errors_("error.dot"),
  write_graphviz_file_on_error_(true),
  evaluation_type_name_(evaluation_type_name),
  sorting_called_(false),
#ifdef PHX_ALLOW_MULTIPLE_EVALUATORS_FOR_SAME_FIELD
  allow_multiple_evaluators_for_same_field_(true),
#else
  allow_multiple_evaluators_for_same_field_(false),
#endif
  build_device_dag_(false),
  field_use_range_evaluated_(false),
  unshared_evaluated_(false)
{ }

//=======================================================================
template<typename Traits>
PHX::DagManager<Traits>::~DagManager()
{ }

//=======================================================================
template<typename Traits>
void PHX::DagManager<Traits>::
requireField(const PHX::FieldTag& t)
{
  FTPredRef pred(t);
  std::vector< Teuchos::RCP<PHX::FieldTag> >::iterator i =
    std::find_if(required_fields_.begin(), required_fields_.end(), pred);

  if (i == required_fields_.end())
    required_fields_.push_back(t.clone());
}

//=======================================================================
template<typename Traits>
void PHX::DagManager<Traits>::
registerEvaluator(const Teuchos::RCP<PHX::Evaluator<Traits> >& p)
{
#ifdef PHX_TEUCHOS_TIME_MONITOR
  // Add counter to name so that all timers have unique names
  static int count=0;
  std::stringstream uniqueName;
  uniqueName << "Phalanx: Evaluator " << count++ <<": [" << evaluation_type_name_ << "] ";
  evalTimers.push_back(
     Teuchos::TimeMonitor::getNewTimer(uniqueName.str() + p->getName()));
#endif

  // insert evaluated fields into map, check for multiple evaluators
  // that provide the same field.
  nodes_.push_back(PHX::DagNode<Traits>(static_cast<int>(nodes_.size()),p));
  const std::vector<Teuchos::RCP<PHX::FieldTag>>& evaluatedFields =
    p->evaluatedFields();
  for (auto i=evaluatedFields.cbegin(); i != evaluatedFields.cend(); ++i) {
    auto check = field_to_node_index_.insert(std::make_pair((*i)->identifier(), static_cast<int>(nodes_.size()-1)));
    if (!allow_multiple_evaluators_for_same_field_) {
      TEUCHOS_TEST_FOR_EXCEPTION(check.second == false,
				 PHX::multiple_evaluator_for_field_exception,
				 *this
				 << "\n\nError: PHX::DagManager::registerEvaluator() - The field \""
				 << (*i)->identifier()
				 << "\" that is evaluated by the evaluator named \""
				 << p->getName()
				 << "\" is already evaluated by another registered evaluator named \""
				 << (nodes_[field_to_node_index_[(*i)->identifier()]]).get()->getName()
				 << "\"."
				 // << " Printing evaluators:\n" << *p << "\n"
				 // << (evaluators_[field_to_index_[(*i)->identifier()]])
				 << std::endl);
    }
  }

  // Insert contributed fields into a separate node list
  const auto& contributedFields = p->contributedFields();
  for (auto i=contributedFields.cbegin(); i != contributedFields.cend(); ++i) {
    auto& set = contributed_field_to_node_index_[(*i)->identifier()];
    set.insert(static_cast<int>(nodes_.size()-1));
  }

}

//=======================================================================
template<typename Traits>
void PHX::DagManager<Traits>::
setDefaultGraphvizFilenameForErrors(const std::string& file_name)
{
  graphviz_filename_for_errors_ = file_name;
}

//=======================================================================
template<typename Traits>
void PHX::DagManager<Traits>::
setWriteGraphvizFileOnError(bool write_file)
{
  write_graphviz_file_on_error_ = write_file;
}

//=======================================================================
template<typename Traits>
void PHX::DagManager<Traits>::
sortAndOrderEvaluators()
{
#ifdef PHX_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Phalanx::SortAndOrderEvaluators"));
#endif

  // *************************
  // Color all nodes white, reset the discovery and final times
  // *************************
  for (auto& n : nodes_)
    n.resetDfsParams(PHX::Color::WHITE);

  topoSortEvalIndex.clear();

  // *************************
  // Insert contributed fields into the field_to_node_index_ if there
  // is no evaluator already assigned. Here we support two cases - (1)
  // there is an evaluator that "evalautes" the field and is already
  // assigned but other evaluators declared as "contributor" for the
  // same field are present and (2) there no evaluators that
  // "evaluate" the field, only "contributors". In the second case, we
  // need to pick one of the contributors and insert it into
  // field_to_node_index_ for the dfs algortihm to work. We also need
  // to remove this from the contributed_field_to_node_index_ so that
  // it does not have a cyclic dependency with itself.
  // *************************
  // for (const auto& contrib_field : contributed_field_to_node_index_) {
  //   const auto& identifier = contrib_field.first;
  //   auto node_index_it = field_to_node_index_.find(identifier);
  //   if (node_index_it == field_to_node_index_.end()) {
  //     field_to_node_index_.insert(std::make_pair(identifier, *contrib_field.second.begin()));
  //     contributed_field_to_node_index_[identifier].erase(field_to_node_index_[identifier]);
  //   }
  // }

  // *************************
  // Loop over required fields
  // *************************
  int time = 0;
  for (const auto& req_field : required_fields_) {

    // Look in evaluated fields first
    auto node_index_it = field_to_node_index_.find(req_field->identifier());
    if (node_index_it != field_to_node_index_.end()) {
      auto& node = nodes_[node_index_it->second];
      if (node.color() == PHX::Color::WHITE)
        dfsVisit(node,time);
    }

    // Add in contributed fields for this tag. There could be multiple
    // contrib or we may not have found any evaluated fields for this
    // tag if they are all contributed.
    auto contrib_field_search = contributed_field_to_node_index_.find(req_field->identifier());
    if (contrib_field_search != contributed_field_to_node_index_.end()) {
      const auto& node_list_to_add = contrib_field_search->second;
      for (auto node_index : node_list_to_add) {
        auto& node = nodes_[node_index];
        if (node.color() == PHX::Color::WHITE)
          dfsVisit(node,time);
      }
    }

    // If no evaluated or contrib found, error out.
    if ( (node_index_it == field_to_node_index_.end()) &&
         (contrib_field_search == contributed_field_to_node_index_.end()) ) {

      if (write_graphviz_file_on_error_)
        this->writeGraphvizFileNew(graphviz_filename_for_errors_, true, true);
      
      TEUCHOS_TEST_FOR_EXCEPTION(node_index_it == field_to_node_index_.end(),
                                 PHX::missing_evaluator_exception,
                                 *this
                                 << "\n\nERROR: The required field \""
                                 << req_field->identifier()
				 << "\" does not have an evaluator. Current "
                                 << "list of Evaluators are printed above this "
                                 << "error message.\n");
    }
  }

  // Create a list of fields to allocate
  fields_.clear();
  for (std::size_t i = 0; i < topoSortEvalIndex.size(); i++) {
    const auto& fields = (nodes_[topoSortEvalIndex[i]]).get()->evaluatedFields();
    fields_.insert(fields_.end(),fields.cbegin(),fields.cend());
  }

  // Contributed fields: If a field is only evaluated with
  // "contributed" fields, then it is not in the list of
  // evaluatedFields and therefore not in the fields to allocate
  // vector above. We need to find such fields and add them to the
  // list of fields to allocate.
  for (std::size_t i = 0; i < topoSortEvalIndex.size(); i++) {
    const auto& contrib_fields = (nodes_[topoSortEvalIndex[i]]).get()->contributedFields();
    for (auto& cfield : contrib_fields) {
      const auto& check_it = std::find_if(fields_.begin(), fields_.end(),[&cfield](const Teuchos::RCP<const PHX::FieldTag>& f)
                                          {
                                            if (f->identifier() == cfield->identifier())
                                              return true;
                                            return false;
                                          }
                                          );
      if (check_it == fields_.end())
        fields_.push_back(cfield);
    }
  }

  this->createEvaluatorBindingFieldMap();

  sorting_called_ = true;
}

//=======================================================================
template<typename Traits>
void PHX::DagManager<Traits>::
dfsVisit(PHX::DagNode<Traits>& node, int& time)
{
  node.setColor(PHX::Color::GREY);
  time += 1;
  node.setDiscoveryTime(time);

  // Add the adjacencies.
  // NOTE: we could do this for all nodes before entering the DFS
  // algorithm, but then the safety check forces users to satisfy
  // dependencies for nodes that may potentially NOT be in the final
  // graph.  So we have to do it here when we know we actually will
  // use the node.
  {
    const auto& req_fields = node.get()->dependentFields();
    for (const auto& field : req_fields) {

      // Look in evalauted fields
      auto node_index_it = field_to_node_index_.find(field->identifier());
      if (node_index_it == field_to_node_index_.end()) {

        // If failed to find, look in contributed fields
        auto contrib_field_search = contributed_field_to_node_index_.find(field->identifier());
        if (contrib_field_search == contributed_field_to_node_index_.end()) {

          if (write_graphviz_file_on_error_)
            this->writeGraphvizFileNew(graphviz_filename_for_errors_, true, true);

          TEUCHOS_TEST_FOR_EXCEPTION(node_index_it == field_to_node_index_.end(),
                                     PHX::missing_evaluator_exception,
                                     *this
                                     << "\n\nERROR: The required field \""
                                     << field->identifier()
                                     << "\" does not have an evaluator. Current "
                                     << "list of Evaluators are printed above this "
                                     << "error message.\n\n"
                                     << "\nPlease inspect the DagManager output above, or \n"
                                     << "visually inspect the error graph that was dumped by \n"
                                     << "running the graphviz dot program on the file\n"
                                     << graphviz_filename_for_errors_ << ": \n\n"
                                     << "dot -Tjpg -o error.jpg "
                                     << graphviz_filename_for_errors_ << "\n\n"
                                     << "The above command generates a jpg file, \"error.jpg\"\n"
                                     << "that you can view in any web browser/graphics program.\n");
        }
      }

      if (node_index_it != field_to_node_index_.end())
        node.addAdjacency(node_index_it->second);

      // For required contributed fields, add the extra evaluators as adjacencies too.
      {
        auto contrib_field_search = contributed_field_to_node_index_.find(field->identifier());
        if (contrib_field_search != contributed_field_to_node_index_.end()) {
          const auto& node_list_to_add = contrib_field_search->second;
          for (auto node_to_add : node_list_to_add)
            node.addAdjacency(node_to_add);
        }
      }
    }

    // For "contributed" fields, if an evaluator exists that also
    // "evalautes" this field, then we assume that the evaluator that
    // "evalautes" the field performs the initialization of the field
    // for the contributions. So the contributed fields must have a
    // dependency on the evalauted field evaluator.
    const auto& contrib_fields = node.get()->contributedFields();
    for (const auto& cfield : contrib_fields) {
      const auto& evaluated_field_search = field_to_node_index_.find(cfield->identifier());
      if (evaluated_field_search != field_to_node_index_.end()) {
        node.addAdjacency(evaluated_field_search->second);
      }
    }

  }

  for (auto& adj_node_index : node.adjacencies()) {

    auto& adj_node = nodes_[adj_node_index];

    if (adj_node.color() == PHX::Color::WHITE) {
      dfsVisit(adj_node,time);
    }
    else if (adj_node.color() == PHX::Color::GREY) {

      std::ostringstream os_adj_node;
      this->printEvaluator(*(adj_node.get()),os_adj_node);
      std::ostringstream os_node;
      this->printEvaluator(*(node.get()),os_node);

      if (write_graphviz_file_on_error_)
	this->writeGraphvizFileNew(graphviz_filename_for_errors_, true, true);

      TEUCHOS_TEST_FOR_EXCEPTION(adj_node.color() == PHX::Color::GREY,
				 PHX::circular_dag_exception,
				 *this
				 << "\n\nERROR: In constructing the directed acyclic graph from \n"
				 << "the node dependencies, a circular dependency has been \n"
				 << "identified. The dependence is injected in going from node:\n\n"
				 << os_adj_node.str()
				 << "\n back to node\n\n"
				 << os_node.str()
				 << "\nPlease inspect the DagManager output above, or \n"
				 << "visually inspect the error graph that was dumped by \n"
				 << "running the graphviz dot program on the file\n"
				 << graphviz_filename_for_errors_ << ": \n\n"
				 << "dot -Tjpg -o error.jpg "
				 << graphviz_filename_for_errors_<< "\n\n"
				 << "The above command generates a jpg file, \"error.jpg\"\n"
				 << "that you can view in any web browser/graphics program.\n");
    }
  }
  node.setColor(PHX::Color::BLACK);
  time += 1;
  node.setFinalTime(time);
  topoSortEvalIndex.push_back(node.index()); // for topo sort
}

//=======================================================================
template<typename Traits>
void PHX::DagManager<Traits>::
printEvaluator(const PHX::Evaluator<Traits>& e, std::ostream& os) const
{
  os << e.getName() << "\n";
  if (e.evaluatedFields().size() > 0) {
    os << "  *Evaluated Fields:\n";
    for (const auto& f : e.evaluatedFields())
      os << "    " << f->identifier() << "\n";
  }
  if (e.contributedFields().size() > 0) {
    os << "  *Contributed Fields:\n";
    for (const auto& f : e.contributedFields())
      os << "    " << f->identifier() << "\n";
  }
  if (e.dependentFields().size() > 0) {
    os << "  *Dependent Fields:\n";
    for (const auto& f : e.dependentFields())
      os << "    " << f->identifier() << "\n";
  }
  if (e.unsharedFields().size() > 0) {
    os << "  *Unshared Fields:\n";
    for (const auto& f : e.unsharedFields())
      os << "    " << f->identifier() << "\n";
  }
}

//=======================================================================
template<typename Traits>
void PHX::DagManager<Traits>::
postRegistrationSetup(typename Traits::SetupData d,
		      PHX::FieldManager<Traits>& vm,
                      const bool& buildDeviceDAG)
{
  // Call each evaluators' post registration setup
  for (std::size_t n = 0; n < topoSortEvalIndex.size(); ++n)
    nodes_[topoSortEvalIndex[n]].getNonConst()->postRegistrationSetup(d,vm);

  build_device_dag_ = buildDeviceDAG;
  if (build_device_dag_) {
    device_evaluators_ = Kokkos::View<PHX::DeviceEvaluatorPtr<Traits>*,PHX::Device>("device_evaluators_",topoSortEvalIndex.size());
    for (std::size_t n = 0; n < topoSortEvalIndex.size(); ++n)
      device_evaluators_(n).ptr = nodes_[topoSortEvalIndex[n]].getNonConst()->createDeviceEvaluator();
  }
}

//=======================================================================
template<typename Traits>
void PHX::DagManager<Traits>::
evaluateFields(typename Traits::EvalData d)
{
  for (std::size_t n = 0; n < topoSortEvalIndex.size(); ++n) {

#ifdef PHX_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor Time(*evalTimers[topoSortEvalIndex[n]]);
#endif

#ifdef PHX_DEBUG
    if (nonnull(start_stop_debug_ostream_)) {
      *start_stop_debug_ostream_ << "Phalanx::DagManager: Starting node: "
                                 << nodes_[topoSortEvalIndex[n]].get()->getName()
                                 << std::endl;
    }

    auto print_fields = std::getenv("PHX_PRINT_FIELDS");
    if (print_fields) {
      std::cout << "************************************************\n"
                << "Printing fields BEFORE executing evaluator: "
                << nodes_[topoSortEvalIndex[n]].get()->getName()
                << "\n  Evaluated Fields:\n";
      auto& eft = nodes_[topoSortEvalIndex[n]].get()->evaluatedFields();
      for (const auto& t : eft)
        std::cout << "    " << t->identifier() << std::endl;
      std::cout << "  Contributed Fields:\n";
      auto& cft = nodes_[topoSortEvalIndex[n]].get()->contributedFields();
      for (const auto& t : cft)
        std::cout << "    " << t->identifier() << std::endl;
      std::cout << "  Dependent Fields:\n";
      auto& dft = nodes_[topoSortEvalIndex[n]].get()->dependentFields();
      for (const auto& t : dft)
        std::cout << "    " << t->identifier() << std::endl;
      std::cout  << "************************************************\n";
      nodes_[topoSortEvalIndex[n]].get()->printFieldValues(std::cout);
    }
#endif

    using clock = std::chrono::steady_clock;
    std::chrono::time_point<clock> start = clock::now();

    typename PHX::Device().fence(); // temporary fence until UVM in evaluateFields fixed

    nodes_[topoSortEvalIndex[n]].getNonConst()->evaluateFields(d);

    nodes_[topoSortEvalIndex[n]].sumIntoExecutionTime(clock::now()-start);

#ifdef PHX_DEBUG
    if (print_fields) {
      std::cout << "************************************************\n"
                << "Printing fields AFTER executing evaluator: "
                << nodes_[topoSortEvalIndex[n]].get()->getName()
                << "\n  Evaluated Fields:\n";
      auto& eft = nodes_[topoSortEvalIndex[n]].get()->evaluatedFields();
      for (const auto& t : eft)
        std::cout << "    " << t->identifier() << std::endl;
      std::cout << "  Contributed Fields:\n";
      auto& cft = nodes_[topoSortEvalIndex[n]].get()->contributedFields();
      for (const auto& t : cft)
        std::cout << "    " << t->identifier() << std::endl;
      std::cout << "  Dependent Fields:\n";
      auto& dft = nodes_[topoSortEvalIndex[n]].get()->dependentFields();
      for (const auto& t : dft)
        std::cout << "    " << t->identifier() << std::endl;
      std::cout  << "************************************************\n";
      nodes_[topoSortEvalIndex[n]].get()->printFieldValues(std::cout);
    }

    if (nonnull(start_stop_debug_ostream_)) {
      *start_stop_debug_ostream_ << "Phalanx::DagManager: Completed node: "
                                 << nodes_[topoSortEvalIndex[n]].getNonConst()->getName()
                                 << std::endl;
    }
#endif

  }
}

//=======================================================================
// Functor for Device DAG support
namespace PHX {

  template<typename Traits>
  struct RunDeviceDag {

    Kokkos::View<PHX::DeviceEvaluatorPtr<Traits>*,PHX::Device> evaluators_;

    // The EvalData may be pass by reference. Remove the reference so
    // that we copy by value to device.
    const typename std::remove_reference<typename Traits::EvalData>::type data_;

    RunDeviceDag(const Kokkos::View<PHX::DeviceEvaluatorPtr<Traits>*,PHX::Device>& evaluators,
                 typename Traits::EvalData data) :
      evaluators_(evaluators),
      data_(data) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) const
    {
      const int num_evaluators = static_cast<int>(evaluators_.extent(0));
      for (int e=0; e < num_evaluators; ++e) {
        evaluators_(e).ptr->prepareForRecompute(team,data_);
        evaluators_(e).ptr->evaluate(team,data_);
	team.team_barrier();
      }
    }

  };

}

//=======================================================================
template<typename Traits>
void PHX::DagManager<Traits>::
evaluateFieldsDeviceDag(const int& work_size,
			const int& team_size,
			const int& vector_size,
                        typename Traits::EvalData d)
{
  TEUCHOS_ASSERT(build_device_dag_);
  //! The parallel_for kernel launch below will not compile on CUDA
  //! unless relocatable device code (RDC) is enabled for the nvcc
  //! compiler. We also want to build and run phalanx without Device
  //! DAG support on CUDA (i.e. RDC off), so this ifdef will hide the
  //! RDC required code.
#if defined(PHX_ENABLE_DEVICE_DAG)
  Kokkos::parallel_for(Kokkos::TeamPolicy<PHX::exec_space>(work_size,team_size,vector_size),
                       PHX::RunDeviceDag<Traits>(device_evaluators_,d));
#else
  TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
    "ERROR: PHX::DagManger::evalauteFieldsDeviceDAG() is experimental and must be enabled at configure time with Phalanx_ENABLE_DEVICE_DAG=ON.");
#endif
}

//=======================================================================
#ifdef PHX_ENABLE_KOKKOS_AMT
template<typename Traits>
void PHX::DagManager<Traits>::
evaluateFieldsTaskParallel(const int& work_size,
			   typename Traits::EvalData d)
{
  using execution_space = PHX::exec_space;
  using memory_space = PHX::Device::memory_space;
  using policy_type = Kokkos::TaskScheduler<execution_space>;

  // Requested the ability to query policy for required sizes of calls
  // to spawn and wait_all. For now hard code to something
  // reasonable!?!
  const unsigned required_memory = 1000000;
  for (std::size_t n = 0; n < topoSortEvalIndex.size(); ++n) {
    // const auto& node = nodes_[topoSortEvalIndex[n]];
    // const auto& adjacencies = node.adjacencies();
  }

  policy_type policy(memory_space(),required_memory);

  // Issue in reusing vector. The assign doesn't like the change of policy.
  //node_futures_.resize(nodes_.size());
  std::vector<Kokkos::Future<void,PHX::exec_space>> node_futures_(nodes_.size());

  for (std::size_t n = 0; n < topoSortEvalIndex.size(); ++n) {

    auto& node = nodes_[topoSortEvalIndex[n]];
    const auto& adjacencies = node.adjacencies();

    // Since this is registered in the order of the topological sort,
    // we know all dependent futures of a node are already
    // constructed.
    std::vector<Kokkos::Future<void,execution_space>> dependent_futures(adjacencies.size());
    auto adj_iterator = adjacencies.cbegin();
    for (std::size_t i=0; i < adjacencies.size(); ++i,++adj_iterator)
      dependent_futures[i] = node_futures_[*adj_iterator];

    auto future = node.getNonConst()->createTask(policy,work_size,dependent_futures,d);
    TEUCHOS_TEST_FOR_EXCEPTION(future.is_null(), std::logic_error,
                               "Error in PHX::DagManager<Traits>::evaluateFieldsTaskParallel():\n"
                               << "The policy is out of memory. Increase the memory pool size!\n");
    node_futures_[topoSortEvalIndex[n]] = future;
  }

  Kokkos::wait(policy);
}
#endif

//=======================================================================
template<typename Traits>
void PHX::DagManager<Traits>::
preEvaluate(typename Traits::PreEvalData d)
{
  for (std::size_t n = 0; n < topoSortEvalIndex.size(); ++n)
    nodes_[topoSortEvalIndex[n]].getNonConst()->preEvaluate(d);
}

//=======================================================================
template<typename Traits>
void PHX::DagManager<Traits>::
postEvaluate(typename Traits::PostEvalData d)
{
  for (std::size_t n = 0; n < topoSortEvalIndex.size(); ++n)
    nodes_[topoSortEvalIndex[n]].getNonConst()->postEvaluate(d);
}

//=======================================================================
template<typename Traits>
void PHX::DagManager<Traits>::
setEvaluationTypeName(const std::string& evaluation_type_name)
{
  evaluation_type_name_ = evaluation_type_name;
}

//=======================================================================
template<typename Traits>
void PHX::DagManager<Traits>::
writeGraphvizFile(const std::string filename,
		  bool writeEvaluatedFields,
		  bool writeDependentFields,
		  bool ) const
{
  writeGraphvizFileNew(filename,writeEvaluatedFields,writeDependentFields);
}

//=======================================================================
template<typename Traits>
void PHX::DagManager<Traits>::
writeGraphvizFileNew(const std::string filename,
		     bool writeEvaluatedFields,
		     bool writeDependentFields) const
{
  std::ofstream ofs;
  ofs.open(filename.c_str());
  writeGraphviz(ofs,writeEvaluatedFields,writeDependentFields);
  ofs.close();
}

//=======================================================================
template<typename Traits>
void PHX::DagManager<Traits>::
writeGraphviz(std::ostream& ofs,
              bool writeEvaluatedFields,
              bool writeDependentFields) const
{
  ofs << "digraph G {\n";

  // This can be called from inside a DFS during an error, so we can't
  // change the DFS node objects when starting a new search. Need to
  // copy the node vector.
  std::vector<PHX::DagNode<Traits>> nodes_copy;
  nodes_copy.insert(nodes_copy.end(),nodes_.cbegin(),nodes_.cend());

  for (auto& n : nodes_copy)
    n.resetDfsParams(PHX::Color::WHITE);

  // Loop over required fields
  int missing_node_index = nodes_.size();
  for (const auto& req_field : required_fields_) {
    // Look in evaluated fields
    auto node_index_it = field_to_node_index_.find(req_field->identifier());
    if (node_index_it != field_to_node_index_.end()) {
      auto& node = nodes_copy[node_index_it->second];
      if (node.color() == PHX::Color::WHITE)
	writeGraphvizDfsVisit(node,
			      nodes_copy,
			      ofs,
			      writeEvaluatedFields,
			      writeDependentFields);
    }

    // Check for contributed fields
    auto contrib_field_search = contributed_field_to_node_index_.find(req_field->identifier());
    if (contrib_field_search != contributed_field_to_node_index_.end()) {
      const auto& node_list_to_add = contrib_field_search->second;
      for (auto node_index : node_list_to_add) {
        auto& node = nodes_copy[node_index];
        if (node.color() == PHX::Color::WHITE)
          writeGraphvizDfsVisit(node,
                                nodes_copy,
                                ofs,
                                writeEvaluatedFields,
                                writeDependentFields);
      }
    }

    if ( (node_index_it == field_to_node_index_.end()) && 
         (contrib_field_search == contributed_field_to_node_index_.end()) ) {
      ofs << missing_node_index
          << " ["  << "fontcolor=\"red\"" << ", label=\"  ** MISSING EVALUATOR **\\n    "
          << req_field->identifier() << "    **** MISSING ****\"]\n";
      missing_node_index += 1;
    }
  }

  ofs << "}";
}

//=======================================================================
template<typename Traits>
void PHX::DagManager<Traits>::
writeGraphvizDfsVisit(PHX::DagNode<Traits>& node,
		      std::vector<PHX::DagNode<Traits>>& nodes_copy,
		      std::ostream& ofs,
		      const bool writeEvaluatedFields,
		      const bool writeDependentFields) const
{
  node.setColor(PHX::Color::GREY);

  // Add valid adjacencies, write node
  {
    std::string font_color = "";
    std::vector<std::string> dependent_field_labels;

    // Dependent field adjacencies
    const auto& req_fields = node.get()->dependentFields();
    for (const auto& field : req_fields) {
      // Look in evaluated fields
      auto node_index_it = field_to_node_index_.find(field->identifier());
      if (node_index_it == field_to_node_index_.end()) {
        // If failed to find, look in contributed fields
        auto contrib_field_search = contributed_field_to_node_index_.find(field->identifier());
        if (contrib_field_search == contributed_field_to_node_index_.end()) {
          font_color = "red";
          std::string dependent_field_label = field->identifier() + "    **** MISSING ****";
          dependent_field_labels.emplace(dependent_field_labels.end(),dependent_field_label);
        }
        else
          dependent_field_labels.push_back(field->identifier());
      }
      else {
        dependent_field_labels.push_back(field->identifier());
        if (node_index_it != field_to_node_index_.end()) {
          node.addAdjacency(node_index_it->second);
        }
      }

      // Add contributed field dependencies
      const auto& contrib_node_index_set_it = contributed_field_to_node_index_.find(field->identifier());
      if (contrib_node_index_set_it != contributed_field_to_node_index_.end()) {
        for (const auto& cnode_index : (*contrib_node_index_set_it).second) {
          node.addAdjacency(cnode_index);
        }
      }
    }

    // For "contributed" fields, if an evaluator exists that also
    // "evalautes" this field, then we assume that the evaluator that
    // "evalautes" the field performs the initialization of the field
    // for the contributions. So the contributed fields must have a
    // dependency on the evalauted field evaluator.
    const auto& contrib_fields = node.get()->contributedFields();
    for (const auto& cfield : contrib_fields) {
      const auto& evaluated_field_search = field_to_node_index_.find(cfield->identifier());
      if (evaluated_field_search != field_to_node_index_.end()) {
        node.addAdjacency(evaluated_field_search->second);
      }
    }

    // Write the node
    ofs << node.index()
	<< " [fontcolor=\"" << font_color
	<< "\", label=\"" << node.get()->getName();
    if (writeEvaluatedFields) {
      ofs << "\\n   Evaluates:";
      const auto& eval_fields = node.get()->evaluatedFields();
      for (const auto& field : eval_fields)
	ofs << "\\n      " << field->identifier();
    }
    if (writeEvaluatedFields) {
      ofs << "\\n   Contributes:";
      for (const auto& field : contrib_fields)
	ofs << "\\n      " << field->identifier();
    }
    if (writeDependentFields) {
      ofs << "\\n   Dependencies:";
      for(const auto& field : dependent_field_labels)
	ofs << "\\n      " << field;
    }
    ofs << "\"]\n";
  }

  // Write edges and trace adjacencies
  for (auto& adj_node_index : node.adjacencies()) {

    auto& adj_node = nodes_copy[adj_node_index];

    if (adj_node.color() == PHX::Color::WHITE) {
      ofs << node.index() << "->" << adj_node.index() << "\n";
      writeGraphvizDfsVisit(adj_node,
			    nodes_copy,
			    ofs,
			    writeEvaluatedFields,
			    writeDependentFields);
    }
    else if (adj_node.color() == PHX::Color::GREY) {
      ofs << node.index() << "->" << adj_node.index() << " [color=red]\n";
    }
    else { // BLACK node
      ofs << node.index() << "->" << adj_node.index() << "\n";
    }
  }

  node.setColor(PHX::Color::BLACK);
}

//=======================================================================
template<typename Traits>
const std::vector<Teuchos::RCP<PHX::FieldTag>>&
PHX::DagManager<Traits>::getFieldTags()
{
  return fields_;
}

//=======================================================================
template<typename Traits>
bool PHX::DagManager<Traits>::sortingCalled() const
{
  return sorting_called_;
}

//=======================================================================
template<typename Traits>
const std::vector<int>&
PHX::DagManager<Traits>::getEvaluatorInternalOrdering() const
{return topoSortEvalIndex;}

//=======================================================================
template<typename Traits>
const std::vector<PHX::DagNode<Traits>>&
PHX::DagManager<Traits>::getDagNodes() const
{
  return nodes_;
}

//=======================================================================
template<typename Traits>
void PHX::DagManager<Traits>::print(std::ostream& os) const
{
  os << "******************************************************" << std::endl;
  os << "PHX::DagManager" << std::endl;
  os << "Evaluation Type = " << evaluation_type_name_ << std::endl;
  os << "******************************************************" << std::endl;

  os << "\n** Starting Required Field List" << std::endl;
  for (std::size_t i = 0; i < required_fields_.size(); i++) {
    os << *(this->required_fields_[i]) << std::endl;
  }
  os << "** Finished Required Field List" << std::endl;

  os << "\n** Starting Registered Field Evaluators" << std::endl;
  for (std::size_t n=0; n < nodes_.size(); ++n) {
    os << "Evaluator[" << n << "]: ";
    this->printEvaluator(*(nodes_[n].get()),os);
  }
  os << "** Finished Registered Field Evaluators" << std::endl;


  os << "\n** Starting Evaluator Order" << std::endl;
  for (std::size_t k = 0; k < topoSortEvalIndex.size(); ++k) {
    os << k << "    " << topoSortEvalIndex[k] << std::endl;
  }
  os << "\nDetails:\n";
  for (std::size_t n = 0; n < topoSortEvalIndex.size(); ++n) {
    os << "Evaluator[" << topoSortEvalIndex[n] << "]: ";
    this->printEvaluator(*(nodes_[topoSortEvalIndex[n]].get()),os);
  }
  os << "** Finished Provider Evaluation Order" << std::endl;

  os << "******************************************************" << std::endl;
  os << "Finished PHX::DagManager" << std::endl;
  os << "Evaluation Type = " << evaluation_type_name_ << std::endl;
  os << "******************************************************" << std::endl;

}

//=======================================================================
template<typename Traits>
void PHX::DagManager<Traits>::
analyzeGraph(double& speedup, double& parallelizability) const
{
  using std::vector;
  using duration = std::chrono::duration<double>;

  duration T_1 = duration(0.0);
  for (std::size_t n = 0; n < topoSortEvalIndex.size(); ++n) {

    const auto& node = nodes_[topoSortEvalIndex[n]];
    T_1 += node.executionTime();

    duration t_data_ready(0.0);
    const auto& adjacencies = node.adjacencies();
    for (const auto& adj_node_index : adjacencies) {
      t_data_ready = std::max(t_data_ready, nodes_[adj_node_index].finishTime());
    }
    const_cast<PHX::DagNode<Traits>&>(node).setStartTime(t_data_ready);
    const_cast<PHX::DagNode<Traits>&>(node).setFinishTime(t_data_ready + node.executionTime());

  }

  duration T_inf = duration(0.0);
  for (std::size_t n = 0; n < topoSortEvalIndex.size(); ++n) {
    const auto& node = nodes_[topoSortEvalIndex[n]];
    T_inf = std::max(T_inf,node.finishTime());
  }

  speedup = T_1.count() / T_inf.count();
  parallelizability = ( 1.0 - (1.0/speedup) )
    / ( 1.0 - (1.0 / static_cast<double>(topoSortEvalIndex.size())) );

}

//=======================================================================
template<typename Traits>
std::vector<Teuchos::RCP<PHX::Evaluator<Traits>>>&
PHX::DagManager<Traits>::getEvaluatorsBindingField(const PHX::FieldTag& ft)
{
  return field_to_evaluators_binding_[ft.identifier()];
}

//=======================================================================
template<typename Traits>
const std::unordered_map<std::string,std::pair<int,int>>&
PHX::DagManager<Traits>::getFieldUseRange()
{
  if (field_use_range_evaluated_)
    return field_use_range_;

  // Initialize the fields just outside valid range for debugging
  for (const auto& f : fields_)
    field_use_range_[f->identifier()] = std::make_pair(topoSortEvalIndex.size(),-1);

  for (int idx=0; idx < static_cast<int>(topoSortEvalIndex.size()); ++idx) {
    const auto& evaluator = nodes_[topoSortEvalIndex[idx]].getNonConst();

    // Evaluated fields
    for (const auto& f : evaluator->evaluatedFields()) {
      auto& range = field_use_range_[f->identifier()];
      range.first = std::min(range.first,idx);
      range.second = std::max(range.first,idx);
    }
    // Contributed fields
    for (const auto& f : evaluator->contributedFields()) {
      auto& range = field_use_range_[f->identifier()];
      range.first = std::min(range.first,idx);
      range.second = std::max(range.first,idx);
    }
    // Dependent fields
    for (const auto& f : evaluator->dependentFields()) {
      auto& range = field_use_range_[f->identifier()];
      range.first = std::min(range.first,idx);
      range.second = std::max(range.first,idx);
    }
  }

  // Required field values must exist to the end of the DAG evalaution
  // since users will access the data after the DAG is finished
  // running.
  for (const auto& f : required_fields_)
    field_use_range_[f->identifier()].second = topoSortEvalIndex.size()-1;

#ifdef PHX_DEBUG
  for (const auto& f : field_use_range_) {
    TEUCHOS_TEST_FOR_EXCEPTION(f.second.first > f.second.first,
                               std::logic_error,
                               "ERROR - PHX::DagManager::getFieldUseRange() - the field "
                               << f.first << " has an invalid use range of [" << f.second.first << ","
                               << f.second.second << "].");
  }
#endif

  field_use_range_evaluated_ = true;
  return field_use_range_;
}
//=======================================================================
template<typename Traits>
const std::unordered_map<std::string,Teuchos::RCP<PHX::FieldTag>>&
PHX::DagManager<Traits>::getUnsharedFields()
{
  if (unshared_evaluated_)
    return unshared_;

  for (int idx=0; idx < static_cast<int>(topoSortEvalIndex.size()); ++idx) {
    const auto& evaluator = nodes_.at(topoSortEvalIndex[idx]).get();
    const auto& e_unshared = evaluator->unsharedFields();
    for (const auto& f : e_unshared)
      unshared_[f->identifier()] = f;
  }

  unshared_evaluated_ = true;
  return unshared_;
}

//=======================================================================
template<typename Traits>
void
PHX::DagManager<Traits>::
printEvaluatorStartStopMessage(const Teuchos::RCP<std::ostream>& ostr)
{
  start_stop_debug_ostream_ = ostr;
}

//=======================================================================
template<typename Traits>
const std::vector<Teuchos::RCP<PHX::FieldTag>>&
PHX::DagManager<Traits>::getRequiredFields() const
{
  return required_fields_;
}
//=======================================================================
template<typename Traits>
const std::unordered_map<std::string,int>&
PHX::DagManager<Traits>::queryRegisteredFields() const
{
  return field_to_node_index_;
}

//=======================================================================
template<typename Traits>
const std::vector<PHX::DagNode<Traits>>&
PHX::DagManager<Traits>::queryRegisteredEvaluators() const
{
  return nodes_;
}

//=======================================================================
template<typename Traits>
void PHX::DagManager<Traits>::createEvaluatorBindingFieldMap()
{
  field_to_evaluators_binding_.clear();

  for (std::size_t n = 0; n < topoSortEvalIndex.size(); ++n) {
    const auto& node = nodes_[topoSortEvalIndex[n]];
    Teuchos::RCP<PHX::Evaluator<Traits>> e = node.getNonConst();

    for (const auto& f : e->evaluatedFields())
      field_to_evaluators_binding_[f->identifier()].push_back(e);

    for (const auto& f : e->contributedFields())
      field_to_evaluators_binding_[f->identifier()].push_back(e);

    for (const auto& f : e->dependentFields())
      field_to_evaluators_binding_[f->identifier()].push_back(e);
  }
}

//=======================================================================
template<typename Traits>
std::ostream&
PHX::operator<<(std::ostream& os, const PHX::DagManager<Traits>& m)
{
  m.print(os);
  return os;
}

//=======================================================================

#endif
