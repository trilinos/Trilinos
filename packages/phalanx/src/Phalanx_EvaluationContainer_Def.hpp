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


#ifndef PHX_EVALUATION_CONTAINER_DEF_HPP
#define PHX_EVALUATION_CONTAINER_DEF_HPP

#include "Teuchos_Assert.hpp"
#include "Phalanx_Traits.hpp"
#include "Phalanx_Evaluator.hpp"
#include "Phalanx_Evaluator_AliasField.hpp"
#include "Phalanx_TypeStrings.hpp"
#include "Phalanx_KokkosViewFactoryFunctor.hpp"
#include <sstream>
#include <stdexcept>

// *************************************************************************
template <typename EvalT, typename Traits>
PHX::EvaluationContainer<EvalT, Traits>::EvaluationContainer() :
  post_registration_setup_called_(false),
  build_device_dag_(false)
{
  this->dag_manager_.setEvaluationTypeName( PHX::typeAsString<EvalT>() );
}

// *************************************************************************
template <typename EvalT, typename Traits> 
PHX::EvaluationContainer<EvalT, Traits>::~EvaluationContainer()
{

}

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
requireField(const PHX::FieldTag& f)
{
  this->dag_manager_.requireField(f);
}

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
registerEvaluator(const Teuchos::RCP<PHX::Evaluator<Traits> >& p)
{
  this->dag_manager_.registerEvaluator(p);
}

// **************************************************************************
template<typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
aliasField(const PHX::FieldTag& aliasedField,
           const PHX::FieldTag& targetField)
{
  // Check that dimensions are the same
  TEUCHOS_TEST_FOR_EXCEPTION(aliasedField.dataLayout() != targetField.dataLayout(),
                             std::runtime_error,
                             "ERROR: The aliased field \"" << aliasedField.identifier()
                             << "\" data layout does not match the target field \""
                             << targetField.identifier() << "\" data layout!\n\n"
                             << "Aliased Field DataLayout:\n" << aliasedField.dataLayout()
                             << "Target Field DataLayout:\n" << targetField.dataLayout() << "\n");
  // Check that the scalar types are the same
  TEUCHOS_TEST_FOR_EXCEPTION(aliasedField.dataTypeInfo() != targetField.dataTypeInfo(),
                             std::runtime_error,
                             "ERROR: The aliased field \"" << aliasedField.identifier()
                             << "\" scalar type does not match the target field \""
                             << targetField.identifier() << "\" scalar type!\n");

  Teuchos::RCP<PHX::Evaluator<Traits>> e =
    Teuchos::rcp(new PHX::AliasField<EvalT,Traits>(aliasedField,targetField));
  this->registerEvaluator(e);
  
  // Store off to assign memory during allocation phase
  aliased_fields_[aliasedField.identifier()] = targetField.identifier();
}

// *************************************************************************
template <typename EvalT, typename Traits> 
void PHX::EvaluationContainer<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData d,
		      PHX::FieldManager<Traits>& fm,
                      const bool& buildDeviceDAG)
{
  // Figure out all evaluator dependencies
  if ( !(this->dag_manager_.sortingCalled()) )
    this->dag_manager_.sortAndOrderEvaluators();


  // Allocate memory for all fields that are needed
  const std::vector< Teuchos::RCP<PHX::FieldTag> >& var_list = 
    this->dag_manager_.getFieldTags();

  std::vector< Teuchos::RCP<PHX::FieldTag> >::const_iterator  var;

  for (var = var_list.begin(); var != var_list.end(); ++var) {
    // skip allocation if this is an aliased field or an unmanaged field
    if ( aliased_fields_.find((*var)->identifier()) == aliased_fields_.end() &&
         unmanaged_fields_.find((*var)->identifier()) == unmanaged_fields_.end() ) {
      typedef typename PHX::eval_scalar_types<EvalT>::type EvalDataTypes;
      Sacado::mpl::for_each_no_kokkos<EvalDataTypes>(PHX::KokkosViewFactoryFunctor<EvalT>(fields_,*(*var),kokkos_extended_data_type_dimensions_));
      
      TEUCHOS_TEST_FOR_EXCEPTION(fields_.find((*var)->identifier()) == fields_.end(),std::runtime_error,
                                 "Error: PHX::EvaluationContainer::postRegistrationSetup(): could not build a Kokkos::View for field named \""
                                 << (*var)->name() << "\" of type \"" << (*var)->dataTypeInfo().name()
                                 << "\" for the evaluation type \"" << PHX::typeAsString<EvalT>() << "\".");
    }
  }

  // Set the unmanaged field memory (not allocated above)
  for (const auto& field : unmanaged_fields_)
    fields_[field.first] = field.second;

  // Assign aliased fields to the target field memory
  for (const auto& field : aliased_fields_)
    fields_[field.first] = fields_[field.second];
  
  // Bind memory to all fields in all required evaluators
  for (const auto& field : var_list)
    this->bindField(*field,fields_[field->identifier()]);

  // This needs to be set before the dag_manager calls each
  // evaluator's postRegistrationSetup() so that the evaluators can
  // use functions that are only valid after post registration setup
  // is called (e.g query for kokkos extended data type dimensions).
  post_registration_setup_called_ = true;

  build_device_dag_ = buildDeviceDAG;  
  // Make sure that device support has been enabled.
#ifndef PHX_ENABLE_DEVICE_DAG
  TEUCHOS_TEST_FOR_EXCEPTION(buildDeviceDAG, std::runtime_error,
                             "ERROR: useDeviceDAG was set to true in call to postRegistrationSetup(), but this feature was not been enabled during configure. Please rebuild with Phalanx_ENABLE_DEVICE_DAG=ON to use this feature.");
#endif
  
  // Allow users to perform special setup. This used to include
  // manually binding memory for all fields in the evaluators via
  // setFieldData(). NOTE: users should not have to bind memory
  // anymore in the postRegistrationSetup() as we now do it for them
  // above.
  this->dag_manager_.postRegistrationSetup(d,fm,build_device_dag_);
}

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
evaluateFields(typename Traits::EvalData d)
{
#ifdef PHX_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( !(this->setupCalled()) , std::logic_error,
		      "You must call postRegistrationSetup() for each evaluation type before calling the evaluateFields() method for that type!");
#endif

  this->dag_manager_.evaluateFields(d);
}

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
evaluateFieldsDeviceDag(const int& work_size,
			const int& team_size,
			const int& vector_size,
			typename Traits::EvalData d)
{
#ifdef PHX_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( !(this->setupCalled()) , std::logic_error,
		      "You must call postRegistrationSetup() for each evaluation type before calling the evaluateFields() method for that type!");
#endif

  this->dag_manager_.evaluateFieldsDeviceDag(work_size,team_size,vector_size,d);
}

// *************************************************************************
#ifdef PHX_ENABLE_KOKKOS_AMT
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
evaluateFieldsTaskParallel(const int& work_size,
			   typename Traits::EvalData d)
{
#ifdef PHX_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( !(this->setupCalled()) , std::logic_error,
		      "You must call post registration setup for each evaluation type before calling the evaluateFields() method for that type!");
#endif

  this->dag_manager_.evaluateFieldsTaskParallel(work_size,d);
}
#endif

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
preEvaluate(typename Traits::PreEvalData d)
{
#ifdef PHX_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( !(this->setupCalled()) , std::logic_error,
		      "You must call post registration setup for each evaluation type before calling the preEvaluate() method for that type!");
#endif

  this->dag_manager_.preEvaluate(d);
}

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
postEvaluate(typename Traits::PostEvalData d)
{
#ifdef PHX_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( !(this->setupCalled()) , std::logic_error,
		      "You must call post registration setup for each evaluation type before calling the postEvaluate() method for that type!");
#endif

  this->dag_manager_.postEvaluate(d);
}

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
setKokkosExtendedDataTypeDimensions(const std::vector<PHX::index_size_type>& dims)
{
  kokkos_extended_data_type_dimensions_ = dims;
}

// *************************************************************************
template <typename EvalT, typename Traits>
const std::vector<PHX::index_size_type> & PHX::EvaluationContainer<EvalT, Traits>::
getKokkosExtendedDataTypeDimensions() const 
{
#ifdef PHX_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( !(this->setupCalled()) , std::logic_error,
		      "You must call post registration setup for each evaluation type before calling the postEvaluate() method for that type!");
#endif
  return kokkos_extended_data_type_dimensions_;
}

// *************************************************************************
template <typename EvalT, typename Traits>
PHX::any
PHX::EvaluationContainer<EvalT, Traits>::getFieldData(const PHX::FieldTag& f)
{
  //return fields_[f.identifier()];
  auto a = fields_.find(f.identifier());
   if (a==fields_.end()){
    std::cout << " PHX::EvaluationContainer<EvalT, Traits>::getFieldData can't find an f.identifier() "<<  f.identifier() << std::endl;
   }
  return a->second;
}

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
setUnmanagedField(const PHX::FieldTag& f, const PHX::any& a)
{
  // An unmanaged field is a MDField where the user has manually
  // allocated the underlying memory for the field. If setup was
  // already called, we need to reassign the memory and rebind all
  // the evalautors that use the unmanaged field. If setup has not
  // been called, then we can store off the memory and assign normally
  // as part of the postRegistrationSetup() process.
  if (this->setupCalled())
    this->bindField(f,a);
  else
    unmanaged_fields_[f.identifier()] = a;

#ifdef PHX_ENABLE_DEVICE_DAG
  if (build_device_dag_) {
    TEUCHOS_ASSERT(false); // need to rebuild device dag
  }
#endif

}

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
bindField(const PHX::FieldTag& f, const PHX::any& a)
{
  auto s = fields_.find(f.identifier());

  if (s == fields_.end()) {
    std::stringstream st;
    st << "\n ERROR in PHX::EvaluationContainer<EvalT, Traits>::bindField():\n" 
       << " Failed to bind field: \"" <<  f.identifier() << "\"\n" 
       << " for evaluation type \"" << PHX::typeAsString<EvalT>() << "\".\n"
       << " This field is not used in the Evaluation DAG.\n";
    
    throw std::runtime_error(st.str());
  }

  // Set the new memory
  fields_[f.identifier()] = a;

  // Loop through evalautors and rebind the field
  auto& evaluators = this->dag_manager_.getEvaluatorsBindingField(f);
  for (auto& e : evaluators)
    e->bindField(f,a);
}

// *************************************************************************
template <typename EvalT, typename Traits>
bool PHX::EvaluationContainer<EvalT, Traits>::setupCalled() const
{
  return post_registration_setup_called_;
}

// *************************************************************************
template <typename EvalT, typename Traits>
const std::string PHX::EvaluationContainer<EvalT, Traits>::
evaluationType() const
{
  return PHX::typeAsString<EvalT>();
}

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::print(std::ostream& os) const
{
  std::string type = PHX::typeAsString<EvalT>();

  os << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  os << "Starting PHX::EvaluationContainer Output" << std::endl;
  os << "Evaluation Type = " << type << std::endl;
  os << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  os << this->dag_manager_ << std::endl;
  for (std::unordered_map<std::string,PHX::any>::const_iterator i = 
	 fields_.begin(); i != fields_.end(); ++i)
    os << i->first << std::endl;
  os << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  os << "Finished PHX::EvaluationContainer Output" << std::endl;
  os << "Evaluation Type = " << type << std::endl;
  os << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  os << std::endl;
}

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
analyzeGraph(double& speedup, double& parallelizability) const
{
  this->dag_manager_.analyzeGraph(speedup,parallelizability);
}

// *************************************************************************
template <typename EvalT, typename Traits>
void
PHX::EvaluationContainer<EvalT, Traits>::buildDag()
{
  this->dag_manager_.sortAndOrderEvaluators();
}

// *************************************************************************
template <typename EvalT, typename Traits>
const std::vector<Teuchos::RCP<PHX::FieldTag>>&
PHX::EvaluationContainer<EvalT, Traits>::getFieldTags()
{
  return this->dag_manager_.getFieldTags();
}

// *************************************************************************
template <typename EvalT, typename Traits>
void
PHX::EvaluationContainer<EvalT, Traits>::
printEvaluatorStartStopMessage(const Teuchos::RCP<std::ostream>& ostr)
{
  this->dag_manager_.printEvaluatorStartStopMessage(ostr);
}

// *************************************************************************

#endif 
