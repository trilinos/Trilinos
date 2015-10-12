// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER


#ifndef PHX_FIELD_EVALUATOR_MANAGER_DEF_HPP
#define PHX_FIELD_EVALUATOR_MANAGER_DEF_HPP

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

#include "boost/version.hpp"
#if defined(BOOST_VERSION)&&(BOOST_VERSION>=104200)
// icpc cannot handle adjacency_list.hpp in BGL after boost 1.56.0
#if (BOOST_VERSION<105600) || !defined(__INTEL_COMPILER)
#include "boost/graph/graphviz.hpp"
#include "boost/tuple/tuple.hpp"
#endif
#endif

//=======================================================================
template<typename Traits>
PHX::EvaluatorManager<Traits>::
EvaluatorManager(const std::string& evaluation_type_name) :
  graphviz_filename_for_errors_("error.dot"),
  write_graphviz_file_on_error_(true),
  evaluation_type_name_(evaluation_type_name),
  sorting_called_(false),
#ifdef PHX_ENABLE_NEW_DFS_ALGORITHM
  use_new_dfs_algorithm_(true),
#else
  use_new_dfs_algorithm_(false),
#endif
#ifdef PHX_ALLOW_MULTIPLE_EVALUATORS_FOR_SAME_FIELD
  allow_multiple_evaluators_for_same_field_(true)
#else
  allow_multiple_evaluators_for_same_field_(false)
#endif
{ }

//=======================================================================
template<typename Traits>
PHX::EvaluatorManager<Traits>::~EvaluatorManager()
{ }

//=======================================================================
template<typename Traits>
void PHX::EvaluatorManager<Traits>::
requireField(const PHX::FieldTag& t)
{
  FTPredRef pred(t);
  std::vector< Teuchos::RCP<PHX::FieldTag> >::iterator i = 
    std::find_if(fields_.begin(), fields_.end(), pred);
  
  if (i == fields_.end()) {
    fields_.push_back(t.clone());
    required_fields_.push_back(t.clone());
  }
}

//=======================================================================
template<typename Traits>
void PHX::EvaluatorManager<Traits>::
registerEvaluator(const Teuchos::RCP<PHX::Evaluator<Traits> >& p)
{
  varProviders.push_back(p);
  providerVariables.push_back(p->evaluatedFields());
  providerRequirements.push_back(p->dependentFields());
  providerNames.push_back(p->getName());

#ifdef PHX_TEUCHOS_TIME_MONITOR
  // Add counter to name so that all timers have unique names
  static int count=0;
  std::stringstream uniqueName;
  uniqueName << "Phalanx: Evaluator " << count++ <<": ";
  evalTimers.push_back(
     Teuchos::TimeMonitor::getNewTimer(uniqueName.str() + p->getName()));
#endif

  // insert evaluated fields into map, check for multiple evaluators
  // that provide the same field.
  nodes_.push_back(PHX::DagNode<Traits>(static_cast<const int>(nodes_.size()),p));
  const std::vector<Teuchos::RCP<PHX::FieldTag>>& evaluatedFields = 
    p->evaluatedFields();
  for (auto i=evaluatedFields.cbegin(); i != evaluatedFields.cend(); ++i) {
    auto check = field_to_node_index_.insert(std::make_pair((*i)->identifier(), static_cast<int>(nodes_.size()-1)));
    if (!allow_multiple_evaluators_for_same_field_) {
      TEUCHOS_TEST_FOR_EXCEPTION(check.second == false,
				 PHX::multiple_evaluator_for_field_exception,
				 *this
				 << "\n\nError: PHX::EvaluatorManager::registerEvaluator() - The field \"" 
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
}

//=======================================================================
template<typename Traits>
void PHX::EvaluatorManager<Traits>::
setDefaultGraphvizFilenameForErrors(const std::string& file_name)
{
  graphviz_filename_for_errors_ = file_name;
}

//=======================================================================
template<typename Traits>
void PHX::EvaluatorManager<Traits>::
setWriteGraphvizFileOnError(bool write_file)
{
  write_graphviz_file_on_error_ = write_file;
}

//=======================================================================
template<typename Traits>
void PHX::EvaluatorManager<Traits>::
sortAndOrderEvaluators()
{
  if (use_new_dfs_algorithm_)
    sortAndOrderEvaluatorsNew();
  else
    sortAndOrderEvaluatorsOld();
}

//=======================================================================
template<typename Traits>
void PHX::EvaluatorManager<Traits>::
sortAndOrderEvaluatorsOld()
{
#ifdef PHX_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Phalanx::SortAndOrderEvaluators"));
#endif

  if (sorting_called_) {
    std::string msg = "Setup was already called.  ";
    msg += "Don't call setup more than once!";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
  }

  // Construct the order in which providers need to be called
  this->createProviderEvaluationOrder();

  /*
   * After we have figured out which providers to call, we need to
   * ensure that all variables a provider will provide are added to
   * the varaible list.  For example, if someone writes a provider
   * that evaluates both DENSITY and VISCOSITY, but the physics only
   * requires DENSITY, an exception would be thrown when the provider
   * tries to get an index for the VISCOSITY variable from the
   * VariableArray.  We need to ensure that all provided variables
   * have storage allocated in the array - i.e. register the VISCOSITY
   * as a variable if it was not.
   */
  for (std::size_t i = 0; i < providerEvalOrderIndex.size(); i++) {
    std::size_t k = providerEvalOrderIndex[i];
    for (std::size_t j = 0; j <providerVariables[k].size(); j++)
      this->requireField(*(providerVariables[k][j]));
  }
  
  sorting_called_ = true;
}

//=======================================================================
template<typename Traits>
void PHX::EvaluatorManager<Traits>::
sortAndOrderEvaluatorsNew()
{
#ifdef PHX_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Phalanx::SortAndOrderEvaluatorsNew"));
#endif
  
  // Color all nodes white, reset the discovery and final times
  for (auto& n : nodes_)
    n.resetDfsParams(PHX::Color::WHITE);

  providerEvalOrderIndex.clear();

  // Loop over required fields
  int time = 0;
  for (const auto& req_field : required_fields_) {

    auto node_index_it = field_to_node_index_.find(req_field->identifier());

    if (node_index_it == field_to_node_index_.end()) {
      
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
    
    auto& node = nodes_[node_index_it->second];
    if (node.color() == PHX::Color::WHITE)
      dfsVisit(node,time);
  }

  // Create a list of fields to allocate
  fields_.clear();
  for (std::size_t i = 0; i < providerEvalOrderIndex.size(); i++) {
    const auto& fields = (nodes_[providerEvalOrderIndex[i]]).get()->evaluatedFields();
    fields_.insert(fields_.end(),fields.cbegin(),fields.cend());
  }

  sorting_called_ = true;
}

//=======================================================================
template<typename Traits>
void PHX::EvaluatorManager<Traits>::
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
      auto node_index_it = field_to_node_index_.find(field->identifier());

      if (node_index_it == field_to_node_index_.end()) {

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
				   << "\nPlease inspect the EvaluatorManager output above, or \n"
				   << "visually inspect the error graph that was dumped by \n"
				   << "running the graphviz dot program on the file\n" 
				   << graphviz_filename_for_errors_ << ": \n\n"
				   << "dot -Tjpg -o error.jpg " 
				   << graphviz_filename_for_errors_ << "\n\n"
				   << "The above command generates a jpg file, \"error.jpg\"\n"
				   << "that you can view in any web browser/graphics program.\n");
      }

      node.addAdjacency(node_index_it->second);
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
				 << "\nPlease inspect the EvaluatorManager output above, or \n"
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
  providerEvalOrderIndex.push_back(node.index()); // for topo sort
}

//=======================================================================
template<typename Traits>
void PHX::EvaluatorManager<Traits>::
printEvaluator(const PHX::Evaluator<Traits>& e, std::ostream& os) const
{
  os << "Name=" << e.getName() << "\n";
  os << "  Evaluated:\n";
  for (const auto& f : e.evaluatedFields()) 
    os << "    " << f->identifier() << "\n";
  os << "  Dependent:\n";
  for (const auto& f : e.dependentFields()) 
    os << "    " << f->identifier() << "\n";
  // os << "Dependent:\n";
  //    <<
}

//=======================================================================
template<typename Traits>
void PHX::EvaluatorManager<Traits>::
postRegistrationSetup(typename Traits::SetupData d,
		      PHX::FieldManager<Traits>& vm)
{
  // Call each providers' post registration setup
  for (std::size_t i = 0; i < providerEvalOrderIndex.size(); i++)
    (varProviders[providerEvalOrderIndex[i]])->postRegistrationSetup(d,vm);
}

//=======================================================================
template<typename Traits>
void PHX::EvaluatorManager<Traits>::createProviderEvaluationOrder()
{
  // Before sorting provider order, we need to add any intermediate
  // variables to the fields_ that are not specified by the operators.
  bool done = false;
  while (!done) {
    bool addedVariables = false;
    
    for (std::size_t i = 0; i < fields_.size(); i++) {
      const PHX::FieldTag& v = *(fields_[i]);
      
      // Loop over providers and add any requirements as variables.
      for (std::size_t prov = 0; prov < providerVariables.size(); prov++) {
	for (std::size_t var = 0; var < providerVariables[prov].size(); var++) {
	  if (*(providerVariables[prov][var]) == v) {
	    // Loop over requirements to see if they are in the variable list.
	    for (std::size_t r = 0; r < providerRequirements[prov].size(); r++) {
	      bool isVariable = false;
	      for (std::size_t j = 0; j < fields_.size(); j++) {
		if (*(fields_[j]) == *(providerRequirements[prov][r]))
		  isVariable = true;
	      }
	      if (!isVariable) {
		fields_.push_back(providerRequirements[prov][r]);
		addedVariables = true;
	      }
	    }
	  }
	}
      }
    }
    if (!addedVariables)
      done = true;
  }
  
  std::vector<Teuchos::RCP<PHX::FieldTag> > tmpList = fields_;
  std::vector<Teuchos::RCP<PHX::FieldTag> > tmpProvided;
  
  // Loop over variable list until it is empty or we fail to remove var
  while (tmpList.size() > 0) {
    
    bool removedVariable = false;
    
    // Loop over all variables still in the list until we find a
    // Provider that can remove a variable
    bool foundProvider = false;
    int providerIndex = -1;
    for (std::size_t var = 0; var < tmpList.size(); var++) {
      
      foundProvider = false;
      providerIndex = -1;
      
      // Loop over variable providers to find one that supplies this variable
      for (std::size_t prov = 0; prov < varProviders.size(); prov++) {
	
	// Loop over provided variable names in provider[prov]
	for (std::size_t i = 0; i < providerVariables[prov].size(); i++) {
	  
	  if (*(tmpList[var]) == *(providerVariables[prov][i])) {
	    foundProvider = true;
	    providerIndex = prov;
	    break;
	  }
	  
	}
	
	if (foundProvider)
	  break;
      }
      

      // Make sure requirements are satisfied for this provider
      bool requirementsSatisfied = true;
      
      if (foundProvider) {
	if (providerRequirements[providerIndex].size() > 0) {
	  
	  for (std::size_t req = 0;
	       req < providerRequirements[providerIndex].size();
	       req++) {
	    bool requiredVariableFound = false;
	    for (std::size_t j = 0; j < tmpProvided.size(); j++) {
	      if (*(providerRequirements[providerIndex][req]) == 
		  *(tmpProvided[j]))
		requiredVariableFound = true;
	    }
	    if (!requiredVariableFound) {
	      requirementsSatisfied = false;
	      break;
	    }
	    
	  }
	}
      }
      
      if (foundProvider && requirementsSatisfied) {
	
	// Remove the variable and exit loop
	std::vector<Teuchos::RCP<PHX::FieldTag> >::iterator p = 
	  tmpList.begin();
	tmpList.erase(p+var);
	// Add all vars to provided list and remove all variables
	// that this provider adds
	for (std::size_t i = 0; i < providerVariables[providerIndex].size(); i++) {
	  tmpProvided.push_back(providerVariables[providerIndex][i]);
	  for (std::size_t j = 0; j < tmpList.size(); j++) {
	    if (*(providerVariables[providerIndex][i]) == *(tmpList[j])) {
	      std::vector<Teuchos::RCP<PHX::FieldTag> >::iterator a = 
		tmpList.begin();
	      tmpList.erase(a+j);
	      break;
	    }
	  }
	}
	providerEvalOrderIndex.push_back(providerIndex);
	removedVariable = true;
	break;
      }

    }  // for (std::size_t var = 0; var < tmpList.size(); var++) {

    if (!removedVariable) {
      std::string msg;

      msg += "\n**************************\n";
      msg += "\nError in EvaluatorManager:\n";
      msg += "\n**************************\n";
      msg += "\nPrinting EvaluatorManager:\n";
      std::ostringstream ost2;
      ost2 << *this << std::endl;
      msg += ost2.str();

      msg += "EvaluatorManager: ";
      msg += evaluation_type_name_;
      msg += " \nCould not meet dependencies!\n";
      msg += "The following variables either have no provider or have a\n";
      msg += "provider but could not satisfy provider requirements:\n\n";
      std::ostringstream ost;
      for (std::size_t i = 0; i < tmpList.size(); i++)
	ost << *(tmpList[i]) << std::endl;
      msg += ost.str();

      msg += "\nPlease look at the EvaluatorManager output above, or \n";
      msg += "visually inspect the error graph that was dumped by \n";
      msg += "running the graphviz dot program on the file error.dot: \n";
      msg += "> dot -Tjpg -o error.jpg error.dot\n\n";
      msg += "The above command generates a jpg file, \"error.jpg\"\n";
      msg += "that you can view in any web browser/graphics program.\n";
	
      std::string filename = "error.dot";
      this->writeGraphvizFile(filename, true, true, true);

      TEUCHOS_TEST_FOR_EXCEPTION(!removedVariable, std::logic_error, msg);
    }
    
  } // While tmpList.size() != 0
  
}

//=======================================================================
template<typename Traits>
void PHX::EvaluatorManager<Traits>::
evaluateFields(typename Traits::EvalData d)
{
  for (std::size_t i = 0; i < providerEvalOrderIndex.size(); i++) {
#ifdef PHX_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor Time(*evalTimers[providerEvalOrderIndex[i]]);
#endif
    (varProviders[providerEvalOrderIndex[i]])->evaluateFields(d);
  }
}

//=======================================================================
template<typename Traits>
void PHX::EvaluatorManager<Traits>::
preEvaluate(typename Traits::PreEvalData d)
{
  for (std::size_t i = 0; i < providerEvalOrderIndex.size(); i++)
    (varProviders[providerEvalOrderIndex[i]])->preEvaluate(d);
}

//=======================================================================
template<typename Traits>
void PHX::EvaluatorManager<Traits>::
postEvaluate(typename Traits::PostEvalData d)
{
  for (std::size_t i = 0; i < providerEvalOrderIndex.size(); i++)
    (varProviders[providerEvalOrderIndex[i]])->postEvaluate(d);
}

//=======================================================================
template<typename Traits>
void PHX::EvaluatorManager<Traits>::
setEvaluationTypeName(const std::string& evaluation_type_name)
{
  evaluation_type_name_ = evaluation_type_name;
}

//=======================================================================
template<typename Traits>
void PHX::EvaluatorManager<Traits>::
writeGraphvizFile(const std::string filename,
		  bool writeEvaluatedFields,
		  bool writeDependentFields,
		  bool debugRegisteredEvaluators) const
{
  
  if (use_new_dfs_algorithm_)
    writeGraphvizFileNew(filename,writeEvaluatedFields,writeDependentFields);
  else
    writeGraphvizFileOld(filename,writeEvaluatedFields,writeDependentFields,debugRegisteredEvaluators);
    
}

//=======================================================================
template<typename Traits>
void PHX::EvaluatorManager<Traits>::
writeGraphvizFileOld(const std::string filename,
		     bool writeEvaluatedFields,
		     bool writeDependentFields,
		     bool debugRegisteredEvaluators) const
{
//#if defined(BOOST_VERSION)&&(BOOST_VERSION>=104200)
// icpc cannot handle adjacency_list.hpp in BGL after boost 1.56.0
#if defined(BOOST_VERSION)&&(BOOST_VERSION>=104200) && ((BOOST_VERSION<105600) || !defined(__INTEL_COMPILER))

  using std::string;
  using std::vector;
  using std::map;
  using std::pair;
  using Teuchos::RCP;
  using PHX::FieldTag;

  TEUCHOS_TEST_FOR_EXCEPTION(!sorting_called_ && !debugRegisteredEvaluators, std::logic_error, "Error sorting of evaluators must be done before writing graphviz file.");

#if (BOOST_VERSION>=104400)
  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
    boost::property<boost::vertex_name_t, std::string, 
    boost::property<boost::vertex_color_t, std::string,
    boost::property<boost::vertex_index_t, std::string> > >,
    boost::property<boost::edge_name_t, std::string> > Graph;
#else
  typedef boost::GraphvizDigraph Graph;
#endif

  Graph g_dot;

#if (BOOST_VERSION>=104400)
  boost::dynamic_properties dp;
  dp.property("label",boost::get(boost::vertex_name, g_dot));
  dp.property("fontcolor",boost::get(boost::vertex_color, g_dot));
  dp.property("id",boost::get(boost::vertex_index, g_dot));
  dp.property("label",boost::get(boost::edge_name, g_dot));
#else
  boost::property_map<Graph,boost::vertex_attribute_t>::type
    vertex_attr_map = get(boost::vertex_attribute, g_dot);
  
  boost::property_map<Graph,boost::edge_attribute_t>::type
    edge_attr_map = get(boost::edge_attribute, g_dot);
#endif

  typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_t;
  typedef typename boost::graph_traits<Graph>::edge_descriptor edge_t;

  // link fields to their evaluators
  std::vector< Teuchos::RCP<PHX::Evaluator<Traits> > > evaluators;
  if (!debugRegisteredEvaluators) {
    for (vector<int>::const_iterator index = providerEvalOrderIndex.begin(); 
	 index != providerEvalOrderIndex.end(); ++index)
      evaluators.push_back(varProviders[*index]);
  }
  else{
    evaluators = varProviders;
  }

  map<string,vertex_t> field_to_evaluator_index;
  vertex_t index = 0;
  for (typename vector< RCP<PHX::Evaluator<Traits> > >::const_iterator 
	 evaluator = evaluators.begin(); evaluator != evaluators.end(); 
       ++evaluator, ++index) {

    const vector< RCP<FieldTag> >& eval_fields = 
      (*evaluator)->evaluatedFields();

    for (vector< RCP<FieldTag> >::const_iterator tag = 
	   eval_fields.begin(); tag != eval_fields.end(); ++tag) {
      
      field_to_evaluator_index[(*tag)->identifier()] = index;

    }
  }

  // Create an edgelist with unique edges (by insterting into a map)
  map<string,pair<vertex_t,vertex_t> > graph_edges;
  for (map<string,std::size_t>::const_iterator field = 
	 field_to_evaluator_index.begin(); 
       field != field_to_evaluator_index.end(); ++field) {

    const vector< RCP<FieldTag> >& dep_fields = 
      (evaluators[field->second])->dependentFields();

    for (vector< RCP<FieldTag> >::const_iterator dep_field = 
	   dep_fields.begin(); dep_field != dep_fields.end(); ++dep_field) {

      // Only add the edge of the out node exists
      map<string,vertex_t>::const_iterator search = 
	field_to_evaluator_index.find((*dep_field)->identifier());
      if (search != field_to_evaluator_index.end()) {

	std::ostringstream edge_name;
	edge_name << field->second << ":" 
		  << field_to_evaluator_index[(*dep_field)->identifier()];
	
	graph_edges[edge_name.str()] = std::pair<vertex_t,vertex_t>
	  (field->second, field_to_evaluator_index[(*dep_field)->identifier()]);
      }
    }
  }


  // Create edge graph between evaluators
  for (map<string,pair<vertex_t,vertex_t> >::const_iterator edge = 
	 graph_edges.begin(); edge != graph_edges.end(); ++edge) {

    std::pair<edge_t, bool> boost_edge = 
      boost::add_edge(edge->second.first, edge->second.second, g_dot);

    
#if (BOOST_VERSION>=104400)
    boost::put("label",dp,boost_edge.first,std::string(edge->first));
#else
    edge_attr_map[boost_edge.first]["label"] = edge->first;
#endif
  }

  boost::graph_traits<Graph>::vertex_iterator vi, vi_end;
  for (boost::tie(vi, vi_end) = vertices(g_dot); vi != vi_end; ++vi) {
    
    string label = evaluators[*vi]->getName();

    if (writeEvaluatedFields) {
    
      const vector< RCP<FieldTag> >& eval_fields = 
	(evaluators[*vi])->evaluatedFields();
      
      label += "\\n   Evaluates:";
      if (eval_fields.size() > 0) {
	for (vector< RCP<FieldTag> >::const_iterator field = 
	       eval_fields.begin(); field != eval_fields.end(); ++field) {
	  label += "\\n     ";
	  label += (*field)->name()  
	    + " : " + (*field)->dataLayout().identifier()
	    + " : " + Teuchos::demangleName((*field)->dataTypeInfo().name());
	}
      }
      else 
	label += " None!";
    
    }

    if (writeDependentFields) {
      const vector< RCP<FieldTag> >& dep_fields = 
	(evaluators[*vi])->dependentFields();

      label += "\\n   Dependencies:";
      if (dep_fields.size() > 0) {
	for (vector< RCP<FieldTag> >::const_iterator field = 
	       dep_fields.begin(); field != dep_fields.end(); ++field) {

	  // Mark any broken evaluators in red
	  bool found = true;
	  if (debugRegisteredEvaluators) {

	    map<string,vertex_t>::const_iterator testing = 
	      field_to_evaluator_index.find((*field)->identifier());
	    if (testing == field_to_evaluator_index.end()) {
	      found = false;
#if (BOOST_VERSION>=104400)
	      boost::put("fontcolor",dp,*vi,std::string("red"));
#else
	      vertex_attr_map[*vi]["fontcolor"] = "red";
#endif
	    }

	  }
	  
	  if (found)
	    label += "\\n     ";
	  else
	    label += "\\n     *****MISSING**** ";
	    
	  label += (*field)->name() 
	    + " : " + (*field)->dataLayout().identifier()
	    + " : " + Teuchos::demangleName((*field)->dataTypeInfo().name());
	}
      }
      else 
	label += " None!";

    }

#if (BOOST_VERSION>=104400)
    boost::put("label",dp,*vi,label);
#else
    vertex_attr_map[*vi]["label"] = label;
#endif
  }

  std::ofstream outfile;
  outfile.open (filename.c_str());
#if (BOOST_VERSION>=104400)
  boost::write_graphviz_dp(outfile, g_dot, dp, std::string("id"));
#else
  boost::write_graphviz(outfile, g_dot);
#endif
  outfile.close();

#else

  std::cout << "WARNING: writeGraphvizFile() was called, but this requires a boost library version 1.42 or higher.  \nPlease rebuild Trilinos with a more recent version of boost." << std::endl; 

#endif // defined(BOOST_VERSION)&&(BOOST_VERSION>=104200)
}

//=======================================================================
template<typename Traits>
void PHX::EvaluatorManager<Traits>::
writeGraphvizFileNew(const std::string filename,
		     bool writeEvaluatedFields,
		     bool writeDependentFields) const
{
  std::ofstream ofs;
  ofs.open(filename.c_str());
  
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
    auto node_index_it = field_to_node_index_.find(req_field->identifier());

    if (node_index_it == field_to_node_index_.end()) {
      ofs << missing_node_index 
	  << " ["  << "fontcolor=\"red\"" << ", label=\"  ** MISSING EVALUATOR **\\n    " 
	  << req_field->identifier() << "    **** MISSING ****\"]\n";     
      missing_node_index += 1;
    }
    else {
      auto& node = nodes_copy[node_index_it->second];
      if (node.color() == PHX::Color::WHITE)
	writeGraphvizDfsVisit(node,
			      nodes_copy,
			      ofs,
			      writeEvaluatedFields,
			      writeDependentFields);
    }
  }

  ofs << "}";
  ofs.close();
}

//=======================================================================
template<typename Traits>
void PHX::EvaluatorManager<Traits>::
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

    const auto& req_fields = node.get()->dependentFields(); 
    for (const auto& field : req_fields) {
      auto node_index_it = field_to_node_index_.find(field->identifier());
      
      // failed to find node
      if (node_index_it == field_to_node_index_.end()) {
	font_color = "red";
	std::string dependent_field_label = field->identifier() + "    **** MISSING ****";
	dependent_field_labels.emplace(dependent_field_labels.end(),dependent_field_label);
      }
      else {
	dependent_field_labels.push_back(field->identifier());	
	node.addAdjacency(node_index_it->second);
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
const std::vector< Teuchos::RCP<PHX::FieldTag> >& 
PHX::EvaluatorManager<Traits>::getFieldTags()
{
  return fields_;
}

//=======================================================================
template<typename Traits>
bool PHX::EvaluatorManager<Traits>::sortingCalled() const
{
  return sorting_called_;
}

//=======================================================================
template<typename Traits>
const std::vector<int>& 
PHX::EvaluatorManager<Traits>::getEvaluatorInternalOrdering() const
{return providerEvalOrderIndex;}

//=======================================================================
template<typename Traits>
void PHX::EvaluatorManager<Traits>::print(std::ostream& os) const
{
  os << "******************************************************" << std::endl;
  os << "PHX::EvaluatorManager" << std::endl;
  os << "Evaluation Type = " << evaluation_type_name_ << std::endl;
  os << "******************************************************" << std::endl;

  os << "\n** Starting Required Field List" << std::endl;
  for (std::size_t i = 0; i < fields_.size(); i++) {
    os << *(this->fields_[i]) << std::endl;
  }
  os << "** Finished Required Field List" << std::endl;

  os << "\n** Starting Registered Field Evaluators" << std::endl;
  for (std::size_t i = 0; i < varProviders.size(); i++) {
    os << "Evaluator[" << i << "]: " << providerNames[i] << std::endl;
    os << "  *Evaluates:" << std::endl;
    for (std::size_t j = 0; j < providerVariables[i].size(); j++)
      os << "    " << *((this->providerVariables[i])[j]) << std::endl;
    os << "  *Dependencies:";
    if (providerRequirements[i].size() == 0) {
      os << " None!" << std::endl;
    }
    else {
      os << std::endl;
      for (std::size_t j = 0; j < providerRequirements[i].size(); j++)
	os << "    " << *((this->providerRequirements[i])[j]) << std::endl;
    }
  }
  os << "** Finished Registered Field Evaluators" << std::endl;


  os << "\n** Starting Evaluator Order" << std::endl;
  for (std::size_t k = 0; k < providerEvalOrderIndex.size(); k++) {
    os << k << "    " << providerEvalOrderIndex[k] << std::endl;
  }
  os << "\nDetails:\n";
  for (std::size_t k = 0; k < providerEvalOrderIndex.size(); k++) {
    int i = providerEvalOrderIndex[k];
    os << "Evaluator[" << i << "]: " << providerNames[i] << std::endl;
    os << "  *Evaluates:" << std::endl;
    for (std::size_t j = 0; j < providerVariables[i].size(); j++)
      os << "    " << *((this->providerVariables[i])[j]) << std::endl;
    os << "  *Dependencies:";
    if (providerRequirements[i].size() == 0) {
      os << " None!" << std::endl;
    }
    else {
      os << std::endl;
      for (std::size_t j = 0; j < providerRequirements[i].size(); j++)
	os << "    " << *((this->providerRequirements[i])[j]) << std::endl;
    }
  }
  os << "** Finished Provider Evaluation Order" << std::endl;

  os << "******************************************************" << std::endl;
  os << "Finished PHX::EvaluatorManager" << std::endl;
  os << "Evaluation Type = " << evaluation_type_name_ << std::endl;
  os << "******************************************************" << std::endl;

}

//=======================================================================
template<typename Traits>
std::ostream&
PHX::operator<<(std::ostream& os, const PHX::EvaluatorManager<Traits>& m)
{
  m.print(os);
  return os;
}

//=======================================================================

#endif
