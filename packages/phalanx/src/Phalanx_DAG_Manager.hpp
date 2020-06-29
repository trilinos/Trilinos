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

#ifndef PHX_DAG_MANAGER_HPP
#define PHX_DAG_MANAGER_HPP

#include <string>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <tuple>
#include <utility> // for std::pair
#include "Teuchos_RCP.hpp"
#include "Phalanx_config.hpp"
#include "Phalanx_FieldTag.hpp"
#include "Phalanx_FieldTag_STL_Functors.hpp"
#include "Phalanx_Evaluator.hpp"
#include "Phalanx_Print.hpp"
#include "Phalanx_DAG_Node.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Kokkos_View.hpp"
#include "Phalanx_DeviceEvaluator.hpp"

#ifdef PHX_ENABLE_KOKKOS_AMT
#include "Kokkos_TaskScheduler.hpp"
#endif

namespace PHX {
  
  template<typename Traits> class FieldManager;

  /*! @brief Class to generate the directed acyclic graph (DAG) for
      evaluation.  Determined which Evaluators should be called and
      the order in which to call them such that all dependencies are
      met with consistency.
   */
  template<typename Traits>
  class DagManager {

  public:

    DagManager(const std::string& evaluator_type_name = "???");
    
    ~DagManager();
    
    //! Require a variable to be evaluated.
    void requireField(const PHX::FieldTag& v);
    
    //! Registers an evaluator with the manager.
    void 
    registerEvaluator(const Teuchos::RCP<PHX::Evaluator<Traits> >& p);
    
    //! Sets the default filename for graphiz file generation for DAG construction errors.
    void setDefaultGraphvizFilenameForErrors(const std::string& file_name);

    //! If set to true, a graphviz file will be written during for DAG construction errors.
    void setWriteGraphvizFileOnError(bool write_file);

    /*! Builds the evaluation DAG.  This should only be called after
      all required fields and evaluators are registered. Must be
      called prior to making calls to postRegistrationSetup(),
      evaluateFields(), preEvaluate(), and postEvaluate().  This can
      be called multiple times to build a new DAG if requirements have
      changed or more evaluators have been added.
    */
    void sortAndOrderEvaluators();
    
    /*! Calls post registration setup on all evaluators.
    */
    void postRegistrationSetup(typename Traits::SetupData d,
			       PHX::FieldManager<Traits>& vm,
                               const bool& buildDeviceDAG);
    
    //! Evaluate the required fields using data parallel evaluation on
    //! topological sort of tasks. Calls parallel_for for each node in
    //! DAG.
    void evaluateFields(typename Traits::EvalData d);

    //! Evaluate the required fields using data parallel evaluation on
    //! topological sort of tasks. Uses Device DAG support, calling a
    //! single parallel_for for the entire DAG. This could be faster
    //! than the call to evaluateFields, but all nodes in the DAG are
    //! restricted to the same work_size. This is intended for CUDA
    //! builds where kernel launch overhead can be significant.
    void evaluateFieldsDeviceDag(const int& work_size,
				 const int& team_size,
				 const int& vector_size,
				 typename Traits::EvalData d);
    
#ifdef PHX_ENABLE_KOKKOS_AMT
    /*! \brief Evaluate the fields using hybrid functional
        (asynchronous multi-tasking) and data parallelism.

      @param work_Size The number of items to divide the parallel work over.
      @param d User defined data.
     */
    void evaluateFieldsTaskParallel(const int& work_size,
				    typename Traits::EvalData d);
#endif

    /*! \brief This routine is called before each residual/Jacobian fill.
      
        This routine is called ONCE on the evaluator before the fill
        loop over elements is started.  This allows us to reset global
        objects between each fill.  An example is to reset an evaluator
        that monitors the maximum grid peclet number in a cell.  This
        call would zero out the maximum for a new fill.
    */
    void preEvaluate(typename Traits::PreEvalData d);
    
    /*! \brief This routine is called after each residual/Jacobian fill.
      
        This routine is called ONCE on the evaluator after the fill
        loop over elements is completed.  This allows us to evaluate
        any post fill data.  An example is to print out some
        statistics such as the maximum grid peclet number in a cell.
    */
    void postEvaluate(typename Traits::PostEvalData d);
    
    void setEvaluationTypeName(const std::string& evaluation_type_name);

    /*! Returns the FieldTags for all fields involved in the
        evaluation. Will return an empty vector unless the user has
        built the DAG using one of the following calls:
        postRegistrationSetup(), postRegistrationSetupForType() or
        buildDagForType().

        WARNING: This is a dangerous power user feature. It returns
        non-const field tags so that the fields can be sized after the
        DAG has been created.
     */
    const std::vector<Teuchos::RCP<PHX::FieldTag>>& getFieldTags();

    /// Returns true if sortAndOrderEvaluators has been called.
    bool sortingCalled() const;

    /// Write the DAG to file in graphviz/dot format. This is the deprecated version.
    void writeGraphvizFile(const std::string filename,
			   bool writeEvaluatedFields,
			   bool writeDependentFields,
			   bool debugRegisteredEvaluators) const;

    /// Write the DAG to file in graphviz/dot format.
    void writeGraphvizFileNew(const std::string filename,
			      bool writeEvaluatedFields,
			      bool writeDependentFields) const;

    /// Write the DAG to std::ostream in graphviz/dot format.
    void writeGraphviz(std::ostream& os,
                       bool writeEvaluatedFields,
                       bool writeDependentFields) const;

    //! Printing
    void print(std::ostream& os) const;

    //! Returns the Topological sort ordering. Used for unit testing.
    const std::vector<int>& getEvaluatorInternalOrdering() const;

    //! Returns the internally registered nodes. Used for unit testing.
    const std::vector<PHX::DagNode<Traits>>& getDagNodes() const;

    /** \brief Returns the speedup and parallelizability of the graph.

	Estimates are based on execution times.  This will return
	garbage if the evaluateFields() call has not been made to log
	execution times.
     */
    void analyzeGraph(double& speedup, double& parallelizability) const;

    /** \brief Returns all evaluators that either evaluate or require
        the given field. This is used to bind memory for unmanaged
        views.

        CAUTION: The returned vector is non-const to rebind memory for
        fields in evaluators. Be careful not to corrupt the actual
        vector.
     */
    std::vector<Teuchos::RCP<PHX::Evaluator<Traits>>>& 
    getEvaluatorsBindingField(const PHX::FieldTag& ft);

    /** \brief Returns the evaluator range that the field needs to exist over.
        
        Once a topological sort of evalautors is performed, we have N
        evalautors in a specific order to traverse for the
        evaluation. Each field is used over a subset of the range of
        evaluators. We can reuse field memory if the use range between
        two fields does not overlap. This function returns the range
        over which each field needs to exist. The MemoryManager will
        use this information when binding fields.

        Function is non-const due to lazy evalaution to construct.

        \returns a map where the key is the field identifier and the
        value is a pair of integers representing the inclusive use
        range [0,N-1] over which the field requires memory.
    */
    const std::unordered_map<std::string,std::pair<int,int>>& getFieldUseRange();

    /** Returns a set of field tags for fields that the user has
        requested to NOT be shared with any other field. Unshared
        fields are used to trade off increased memory use for a
        reduction in flops for an evalautor. Unshared fields are a
        corner case where the user can leverage special knowledge
        about how data in a field changes across evaluations. One
        example use case is for FAD types during a Gather operation,
        where we know the off diagonal entries are always zero. The
        evaluator can zero out the FAD array during initialization and
        only change the diagonal (seed value) during an evalaution.
     */
    const std::unordered_map<std::string,Teuchos::RCP<PHX::FieldTag>>& getUnsharedFields();

    /** \brief Print to user specified ostream when each evaluator
        starts and stops. Useful for debugging. Enabled only in debug
        builds.

        @param [in] ostr RCP to output stream. If set to null, this disables printing.
    */
    void printEvaluatorStartStopMessage(const Teuchos::RCP<std::ostream>& ostr);

    /// Returns all fields that the user requested to to be evaluated by the field manager.
    const std::vector<Teuchos::RCP<PHX::FieldTag>>& getRequiredFields() const;

  protected:

    /*! @brief Depth-first search algorithm. */ 
    void dfsVisit(PHX::DagNode<Traits>& node, int& time);
        
    /*! @brief Depth-first search algorithm specialized for writing graphviz output. */ 
    void writeGraphvizDfsVisit(PHX::DagNode<Traits>& node,
			       std::vector<PHX::DagNode<Traits>>& nodes_copy,
			       std::ostream& os,
			       const bool writeEvaluatedFields,
			       const bool writeDependentFields) const;
        
    //! Helper function.
    void printEvaluator(const PHX::Evaluator<Traits>& e, std::ostream& os) const;

    void createEvaluatorBindingFieldMap();
    
  protected:

    //! Fields required by the user.
    std::vector<Teuchos::RCP<PHX::FieldTag>> required_fields_;

    /*! @brief Vector of all registered evaluators. 

      This list may include more nodes than what is needed for the DAG
      evaluation of required fields.
    */
    std::vector<PHX::DagNode<Traits>> nodes_;

    //! Hash map of field key to evaluator index.
    std::unordered_map<std::string,int> field_to_node_index_;

    //! Hash map of contributed field key to evaluator index.
    std::unordered_map<std::string,std::unordered_set<int>>
      contributed_field_to_node_index_;
    
    //! All fields that are needed for the evaluation.
    std::vector< Teuchos::RCP<PHX::FieldTag> > fields_;

    // Timers used when configured with Phalanx_ENABLE_TEUCHOS_TIME_MONITOR.
    std::vector<Teuchos::RCP<Teuchos::Time> > evalTimers;

    /*! @name Evaluation Order Objects
      
        Stores results from a topological sort on the evaluator DAG:
        the order to call evaluators to evaluate fields correctly.
    */
    std::vector<int> topoSortEvalIndex;

    //! Use this name for graphviz file output for DAG construction errors.
    std::string graphviz_filename_for_errors_;

    //! If set to true, will write graphviz file for DAG construction errors.
    bool write_graphviz_file_on_error_;

    std::string evaluation_type_name_;

    //! Flag to tell the setup has been called.
    bool sorting_called_;

    //! Backwards compatibility option: set to true to disable a check that throws if multiple registered evaluators can evaluate the same field. Original DFS algortihm allowed this.  Refactor checks and throws.   
    bool allow_multiple_evaluators_for_same_field_;

#ifdef PHX_ENABLE_KOKKOS_AMT
    //std::vector<Kokkos::Experimental::Future<void,PHX::exec_space>> node_futures_;
#endif

    //! A map that returns all evaluators that bind the memory of a particular field. Key is unique field identifier.  
    std::unordered_map<std::string,std::vector<Teuchos::RCP<PHX::Evaluator<Traits>>>> field_to_evaluators_binding_;

    //! If set to true, allocated DeviceEvaluators for Device DAG for evaluation
    bool build_device_dag_;
    
    //! Contians pointers to DeviceEvaluators for Device DAG support.
    Kokkos::View<PHX::DeviceEvaluatorPtr<Traits>*,PHX::Device> device_evaluators_;

    //! If non-null, in debug builds, the DAG manager will print when an evaluator starts and stops.
    Teuchos::RCP<std::ostream> start_stop_debug_ostream_;

    //! Field use range for topologically sorted evalautors. Key is field identifier, value is inclusive start/stop range.
    std::unordered_map<std::string,std::pair<int,int>> field_use_range_;

    //! True if the field use range has been evaluated.
    bool field_use_range_evaluated_;

    //! Fields the user has requested to NOT share memory.
    std::unordered_map<std::string,Teuchos::RCP<PHX::FieldTag>> unshared_;

    //! True if the unshared fields have been evaluated.
    bool unshared_evaluated_;
  };
  
  template<typename Traits>
  std::ostream& operator<<(std::ostream& os, 
			   const PHX::DagManager<Traits>& m);

}

#include "Phalanx_DAG_Manager_Def.hpp"

#endif
