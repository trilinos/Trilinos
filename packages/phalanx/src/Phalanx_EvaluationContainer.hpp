// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_SCALAR_CONTAINER_HPP
#define PHX_SCALAR_CONTAINER_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Phalanx_config.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_EvaluationContainer_Base.hpp"
#include "Phalanx_FieldTag.hpp"
#include "Phalanx_Evaluator.hpp"
#include <any>
#include <unordered_map>
#include <string>
#include <memory>

namespace PHX {

  class MemoryManager;

  /*! \brief Container that holds all data associated with an evaluation type.

    Handles allocation and binding of all field memory.
  */
  template <typename EvalT, typename Traits>
  class EvaluationContainer : public PHX::EvaluationContainerBase<Traits> {

  public:

    EvaluationContainer();

    ~EvaluationContainer();

    //! Requests that the container must compute this field.
    void requireField(const PHX::FieldTag& f);

    void aliasField(const PHX::FieldTag& aliasedField,
                    const PHX::FieldTag& targetField);

    void
    registerEvaluator(const Teuchos::RCP<PHX::Evaluator<Traits> >& p);

    std::any getFieldData(const PHX::FieldTag& f);

    /** \brief Set the memory for an unmanaged field
     *
     * NOTE: If this method is called after postRegistrationSetup(),
     * the field might be reported as shared when priting even though
     * it is no longer shared (now points to user supplied
     * memory). Output from DAG may be incorrect. Searching the field
     * lists for potential sharing wastes time as this function may be
     * called in the middle of an evaluation, so we will not clean up
     * output or add this to the unmanaged field list unless the user
     * explicitly asks for this cleanup to happen. Execution will
     * always be correct.
    */
    void setUnmanagedField(const PHX::FieldTag& f,
                           const std::any& a,
                           const bool cleanup_output = true);

    //! Bind the memory pointer for a field in all evaluators
    void bindField(const PHX::FieldTag& f, const std::any& a);

    void postRegistrationSetup(typename Traits::SetupData d,
			       PHX::FieldManager<Traits>& fm,
                               const bool& buildDeviceDAG,
                               const bool& minimizeDAGMemoryUse,
                               const PHX::MemoryManager* const memoryManager);

    void evaluateFields(typename Traits::EvalData d);

    void evaluateFieldsDeviceDag(const int& work_size,
				 const int& team_size,
				 const int& vector_size,
				 typename Traits::EvalData d);

#ifdef PHX_ENABLE_KOKKOS_AMT
    /*! \brief Evaluate the fields using hybrid functional (asynchronous multi-tasking) and data parallelism.

      @param work_size The number of work units to parallelize over.
      @param d User defined data.
     */
    void evaluateFieldsTaskParallel(const int& work_size,
				    typename Traits::EvalData d);
#endif

    void preEvaluate(typename Traits::PreEvalData d);

    void postEvaluate(typename Traits::PostEvalData d);

    void setKokkosExtendedDataTypeDimensions(const std::vector<PHX::index_size_type>& dims);

    const std::vector<PHX::index_size_type> & getKokkosExtendedDataTypeDimensions() const;

    //! Return true if the postRegistrationSetupMethod has been called
    bool setupCalled() const;

    const std::string evaluationType() const;

    void print(std::ostream& os) const;

    void analyzeGraph(double& speedup, double& parallelizability) const;

    /*! Build the DAG. This is automatically called by the
        postRegistrationSetup() method. This function is a power user
        feature that allows for cases where the user would like to
        build the dag and query it to use information from the DAG
        prior to allocating and binding the memory to fields.
     */
    void buildDag();

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

    /** \brief Print to user specified ostream when each evaluator
        starts and stops. Useful for debugging. Enabled only in debug
        builds.

        @param [in] ostr RCP to output stream. If set to null, this disables printing.
    */
    void printEvaluatorStartStopMessage(const Teuchos::RCP<std::ostream>& ostr);

    /** Returns the underlying DAGManager. Used for queries, debugging and unit testing. */
    const PHX::DagManager<Traits>& getDagManager() const;

  protected:

    void assignSharedFields();

    bool post_registration_setup_called_;

    std::unordered_map<std::string,std::any> fields_;

    std::unordered_map<std::string,std::any> unmanaged_fields_;

    std::unordered_map<std::string,std::string> aliased_fields_;

    /** Shared fields are fields where their use range in the
        topological sort of the dag does not overlap. Therefore, the
        fields can share the same memory allocation tracker. The key
        is the identifier for the field that will not be allocated
        since it will use another field's memory. The value is a pair
        where first is an RCP to the shared field tag, and second is
        the field string identifier whose memory the shared field will
        point to.
     */
    std::unordered_map<std::string,std::pair<Teuchos::RCP<PHX::FieldTag>,std::string>> shared_fields_;

    std::vector<PHX::index_size_type> kokkos_extended_data_type_dimensions_;

    bool build_device_dag_;

    // Enables shared memory use if set to true.
    bool minimize_dag_memory_use_;

    std::shared_ptr<PHX::MemoryManager> memory_manager_;

    /// Size in bytes of view allocation. This includes padding if the view supports/requires it.
    std::unordered_map<std::string,std::size_t> field_allocation_sizes_;

    std::vector<std::pair<std::size_t,Teuchos::RCP<PHX::FieldTag>>> fields_to_allocate_;
  };

}

#include "Phalanx_EvaluationContainer_Def.hpp"

#endif
