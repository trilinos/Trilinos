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


#ifndef PHX_SCALAR_CONTAINER_HPP
#define PHX_SCALAR_CONTAINER_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Phalanx_config.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_EvaluationContainer_Base.hpp"
#include "Phalanx_FieldTag.hpp"
#include "Phalanx_Evaluator.hpp"
#include "Phalanx_any.hpp"
#include <unordered_map>
#include <string>

namespace PHX {

  /*! \brief Container that holds all data associated with an evaluation type.


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

    PHX::any getFieldData(const PHX::FieldTag& f);

    //! Set the memory for an unmanaged field
    void setUnmanagedField(const PHX::FieldTag& f, const PHX::any& a);

    //! Bind the memory pointer for a field in all evaluators
    void bindField(const PHX::FieldTag& f, const PHX::any& a);

    void postRegistrationSetup(typename Traits::SetupData d,
			       PHX::FieldManager<Traits>& fm,
                               const bool& buildDeviceDAG);

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

    /** \brief Print to user specified ostream when each evalautor
        starts and stops. Useful for debugging. Enabled only in debug
        builds.

        @param [in] ostr RCP to output stream. If set to null, this disables printing.
    */
    void printEvaluatorStartStopMessage(const Teuchos::RCP<std::ostream>& ostr);

  protected:

    bool post_registration_setup_called_;

    std::unordered_map<std::string,PHX::any> fields_;

    std::unordered_map<std::string,PHX::any> unmanaged_fields_;

    std::unordered_map<std::string,std::string> aliased_fields_;
    
    std::vector<PHX::index_size_type> kokkos_extended_data_type_dimensions_;

    bool build_device_dag_;
  };
  
} 

#include "Phalanx_EvaluationContainer_Def.hpp"

#endif 
