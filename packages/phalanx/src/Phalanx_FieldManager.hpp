// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_FIELD_MANAGER_HPP
#define PHX_FIELD_MANAGER_HPP

#include <cstddef>
#include <string>
#include <map>
#include <iostream>
#include <vector>
#include <algorithm>
#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_FieldTag.hpp"
#include "Phalanx_EvaluationContainer_TemplateManager.hpp"

// *******************************
// Forward declarations
// *******************************
namespace PHX {
  template<typename Scalar,typename...Props> class MDField;
  template<typename DataT,int Rank,typename Layout> class Field;
}
namespace Kokkos {
  template<typename DataT,typename... Props> class View;
}

// *******************************
// Field Manager
// *******************************
namespace PHX {

  template<typename Traits>
  class FieldManager {

  public:

    typedef typename PHX::EvaluationContainer_TemplateManager<Traits>::iterator iterator;

    FieldManager();

    ~FieldManager();

    void requireFieldForAllEvaluationTypes(const PHX::FieldTag& t);

    template<typename EvalT>
    void requireField(const PHX::FieldTag& t);

    void registerEvaluatorForAllEvaluationTypes(const Teuchos::RCP< PHX::Evaluator<Traits> >& e);

    template<typename EvalT>
    void registerEvaluator(const Teuchos::RCP< PHX::Evaluator<Traits> >& e);

    void registerEvaluator(typename PHX::FieldManager<Traits>::iterator it,
			   const Teuchos::RCP< PHX::Evaluator<Traits> >& e);

    template<typename EvalT, typename DataT,typename...Props>
    void getFieldData(PHX::MDField<DataT,Props...>& f);

    template<typename EvalT, typename DataT,typename...Props>
    void getFieldData(PHX::MDField<const DataT,Props...>& f);

    template<typename EvalT, typename DataT, int Rank, typename Layout>
    void getFieldData(PHX::Field<DataT,Rank,Layout>& f);

    template<typename EvalT, typename DataT, int Rank, typename Layout>
    void getFieldData(PHX::Field<const DataT,Rank,Layout>& f);

    template<typename EvalT, typename DataT,typename Layout>
    void getFieldData(const PHX::FieldTag& ft,
                      Kokkos::View<DataT,Layout,PHX::Device>& f);

    /*! \brief Allows the user to manage the memory allocation of a
        particular field and dynamically set/swap the memory at any
        time.

        This overrides the field allocated to this array in the
        FieldManager. The fieldManager then sets this new memory
        pointer in all evaluator fields that use it.

        NOTE: this is a very dangerous power user capability as the
        user must allocate the field correctly (remember Sacado AD
        types must have the extra dimensions sized correctly).

        \param cleanup_output (bool) This flag only matters if this
        function is called after postRegistrationSetup() is called. If
        set to true and called after postRegistrationSetup(), this
        will take more execution time to search field lists to cleanup
        data structures for output information. This is important
        because a user could toggle a field that was tagged as shared
        during postRegistrationSetup() into an unmanaged state. The
        code will always perform correctly, but output from this
        object might be confusing as is could report an unmanaged
        field as being shared. We allow users set this flag to false
        and to leave the output in a bad state since they might want
        to call this many times in the middle of an evaluation.
    */
    template<typename EvalT, typename DataT, typename...Props>
    void setUnmanagedField(PHX::MDField<DataT,Props...>& f,
                           const bool cleanup_output = true);

    /*! \brief Allows the user to manage the memory allocation of a
        particular field and dynamically set/swap the memory at any
        time.

        This overrides the field allocated to this array in the
        FieldManager. The fieldManager then sets this new memory
        pointer in all evaluator fields that use it.

        NOTE: this is a very dangerous power user capability as the
        user must allocate the field correctly (remember Sacado AD
        types must have the extra dimensions sized correctly).

        \param cleanup_output (bool) This flag only matters if this
        function is called after postRegistrationSetup() is called. If
        set to true and called after postRegistrationSetup(), this
        will take more execution time to search field lists to cleanup
        data structures for output information. This is important
        because a user could toggle a field that was tagged as shared
        during postRegistrationSetup() into an unmanaged state. The
        code will always perform correctly, but output from this
        object might be confusing as is could report an unmanaged
        field as being shared. We allow users set this flag to false
        and to leave the output in a bad state since they might want
        to call this many times in the middle of an evaluation.
    */
    template<typename EvalT, typename DataT, int Rank, typename Layout>
    void setUnmanagedField(PHX::Field<DataT,Rank,Layout>& f,
                           const bool cleanup_output = true);

    /*! \brief Allows the user to manage the memory allocation of a
        particular field and dynamically set/swap the memory at any
        time.

        This overrides the field allocated to this array in the
        FieldManager. The fieldManager then sets this new memory
        pointer in all evaluator fields that use it.

        NOTE: this is a very dangerous power user capability as the
        user must allocate the field correctly (remember Sacado AD
        types must have the extra dimensions sized correctly).
        \param cleanup_output (bool) This flag only matters if this
        function is called after postRegistrationSetup() is called. If
        set to true and called after postRegistrationSetup(), this
        will take more execution time to search field lists to cleanup
        data structures for output information. This is important
        because a user could toggle a field that was tagged as shared
        during postRegistrationSetup() into an unmanaged state. The
        code will always perform correctly, but output from this
        object might be confusing as is could report an unmanaged
        field as being shared. We allow users set this flag to false
        and to leave the output in a bad state since they might want
        to call this many times in the middle of an evaluation.
    */
    template<typename EvalT, typename DataT, typename Layout>
    void setUnmanagedField(const FieldTag& ft,
                           Kokkos::View<DataT,Layout,PHX::Device>& f,
                           const bool cleanup_ouput = true);

    /*! \brief Makes two fields point to (alias) the same memory for all evaluation types.

       WARNING: this is a very dangerous power user capability. This
       allows users to tell the FieldManager to create a new field
       that points to the same underlying memory as another field. The
       user must be sure that the DataLayouts and Scalar types are the
       same. Only use this BEFORE postRegistrationSetup() is
       called. This injects extra dependencies that must be accounted
       for during DAG construction.

       This is intended for the use case where a user wants to reuse
       an evaluator with hard coded field names but would like to
       rename the evaluated fields without adding naming logic to the
       evaluator.

       @param aliasedField Field that is aliased to the target field's memory
       @param targetField Field whos memory is pointed to by the aliased field
     */
    void aliasFieldForAllEvaluationTypes(const PHX::FieldTag& aliasedField,
                                         const PHX::FieldTag& targetField);

    /*! \brief Makes two fields point to (alias) the same memory for a specific evaluation type.

       WARNING: this is a very dangerous power user capability. This
       allows users to tell the FieldManager to create a new field
       that points to the same underlying memory as another field. The
       user must be sure that the DataLayouts and Scalar types are the
       same. Only use this BEFORE postRegistrationSetup() is
       called. This injects extra dependencies that must be accounted
       for during DAG construction.

       This is intended for the use case where a user wants to reuse
       an evaluator with hard coded field names but would like to
       rename the evaluated fields without adding naming logic to the
       evaluator.

       @param aliasedField Field that is aliased to the target field's memory
       @param targetField Field whos memory is pointed to by the aliased field
     */
    template<typename EvalT>
    void aliasField(const PHX::FieldTag& aliasedField,
                    const PHX::FieldTag& targetField);

    /*! \brief Builds DAG (if not already built) and allocates memory for a single evaluation type

       @param[in] d User defined setup data.
       @param[in] buildDeviceDAG (optional) If set to true, the dag is built on device.
       @param[in] minimizeDAGMemoryUse (optional) If set to true, field memory will be reused in a DAG by binding the same kokkos allocation trackers to non-overlapping fields when possible.
       @param[in] memoryManager (optional) If non-null, field memory allocations will use the memoryManager. This can allow multiple DAGs within a FieldManager and multiple FieldManagers to share/reuse field memory.
     */
    template<typename EvalT>
    void postRegistrationSetupForType(typename Traits::SetupData d,
                                      const bool& buildDeviceDAG = false,
                                      const bool& minimizeDAGMemoryUse = false,
                                      const PHX::MemoryManager* const memoryManager = nullptr);

    /*! \brief Builds DAG (if not already built) and allocates memory for all evaluation types

       @param[in] d User defined setup data.
       @param[in] buildDeviceDAG (optional) If set to true, the dag is built on device.
       @param[in] minimizeDAGMemoryUse (optional) If set to true, field memory will be reused in a DAG by binding the same kokkos allocation trackers to non-overlapping fields when possible.
       @param[in] memoryManager (optional) If non-null, field memory allocations will use the memoryManager. This can allow multiple DAGs within a FieldManager and multiple FieldManagers to share/reuse field memory.
     */
    void postRegistrationSetup(typename Traits::SetupData d,
                               const bool& buildDeviceDAG = false,
                               const bool& minimizeDAGMemoryUse = false,
                               const PHX::MemoryManager* const memoryManager = nullptr);

    //! Evalaute fields with a separate parallel_for for each node in the DAG.
    template<typename EvalT>
    void evaluateFields(typename Traits::EvalData d);

    //! Evalaute fields using Device DAG capability where a single parallel_for evaluates the entire DAG.
    template<typename EvalT>
    void evaluateFieldsDeviceDag(const int& work_size, const int& team_size, const int& vector_size, typename Traits::EvalData d);

#ifdef PHX_ENABLE_KOKKOS_AMT
    /*! \brief Evaluate the fields using hybrid functional (asynchronous multi-tasking) and data parallelism.

      @param work_size The number of parallel work units.
      @param d User defined data.
     */
    template<typename EvalT>
    void evaluateFieldsTaskParallel(const int& work_size,
				    typename Traits::EvalData d);
#endif

    template<typename EvalT>
    void preEvaluate(typename Traits::PreEvalData d);

    template<typename EvalT>
    void postEvaluate(typename Traits::PostEvalData d);

    template<typename EvalT>
    void setKokkosExtendedDataTypeDimensions(const std::vector<PHX::index_size_type>& dims);

    template<typename EvalT>
    const std::vector<PHX::index_size_type>& getKokkosExtendedDataTypeDimensions() const;

    //! Returns DagManager for an evaluation type. Used for query, debug and unit testing.
    template<typename EvalT>
    const PHX::DagManager<Traits>& getDagManager() const;

    //! Return iterator to first EvaluationContainer
    typename FieldManager::iterator begin();

    //! Return iterator to last EvaluationContainer
    typename FieldManager::iterator end();

    //! Writes graphviz dot file for the evaluation type
    template<typename EvalT>
    void writeGraphvizFile(const std::string filename = "graph.dot",
			   bool writeEvaluatedFields = true,
			   bool writeDependentFields = true,
			   bool debugRegisteredEvaluators = false) const;

    //! Writes graphviz dot file for all evaluation types (adds eval type to filename).
    void writeGraphvizFile(const std::string base_filename = "graph",
			   const std::string file_extension = ".dot",
			   bool writeEvaluatedFields = true,
			   bool writeDependentFields = true,
			   bool debugRegisteredEvaluators = false) const;

    void print(std::ostream& os) const;

    template<typename EvalT>
    void analyzeGraph(double& speedup, double& parallelizability) const;

    /*! Builds the DAG for the evalaution type. This should only be
        called after all evaluators are registered and all required
        fields are requested. This method is for power users
        only. This is automatically called during
        postRegistrationSetup() and normally does not have to be
        called by the users. This method allows users to build the DAG
        but then perform other activities prior to allocating the
        fields. An example use case is to delay the sizing of the
        fields in the DataLayouts until right before allocation. The
        user could create the dag and access a list of required fields
        and then do sizing based on information aboutrequired fields.
    */
    template<typename EvalT>
    void buildDagForType();

    /*! Returns the FieldTags for all fields involved in the
        evaluation. Will return an empty vector unless the user has
        built the DAG using one of the following calls:
        postRegistrationSetup(), postRegistrationSetupForType() or
        buildDagForType().

        WARNING: This is a dangerous power user feature. It returns
        non-const field tags so that the fields can be sized after the
        DAG has been created.
     */
    template<typename EvalT>
    const std::vector<Teuchos::RCP<PHX::FieldTag>>&
    getFieldTagsForSizing();

    /** \brief Print to user specified ostream when each evaluator
        starts and stops. Useful for debugging. Enabled only in debug
        builds.

        @param [in] ostr RCP to output stream. If set to null, this disables printing.
    */
    template<typename EvalT>
    void printEvaluatorStartStopMessage(const Teuchos::RCP<std::ostream>& ostr);

  private:

    using SCTM = PHX::EvaluationContainer_TemplateManager<Traits>;

    std::size_t m_num_evaluation_types;

    PHX::EvaluationContainer_TemplateManager<Traits> m_eval_containers;

  };

  template<typename Traits>
  std::ostream& operator<<(std::ostream& os,
			   const PHX::FieldManager<Traits>& vm);

}

#include "Phalanx_FieldManager_Def.hpp"

#endif
